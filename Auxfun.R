#Fast column multiplication
cprod <- function(Y){
	prod <- Y[,1]
	for (j in 2:ncol(Y)){
		prod <- Y[,j]*prod
	}
	return(prod)
	}
#Legendre polynomial
legendre <- function(n, x){
	stopifnot(is.numeric(x), is.numeric(n),
	          length(n) == 1, floor(n) == ceiling(n), n >= 0)
	#Normalize x to [-1,1]
	x <- 2*(x-min(x))/(max(x)-min(x))-1
	x <- as.matrix(x)
	N <- nrow(x)
	if (n == 0) return(rep(1, N))
    # generate the Legendre polynomials up to degree n
    Lp <- matrix(0, nrow=nrow(x), ncol=n+1)
    Lp[,1] <- 1; Lp[,2] <- x
    if (n>1){
        for (i in 3:(n+1)) {
	        j <- i-1
	        k <- i-2
	        Lp[,i] <- (2*k+1)*x*Lp[,j]/j-k*Lp[,k]/j
        }
    }
    lp <- Lp[,n+1]
    return(lp)
}
#Hermite polynomial
hermite <- function(n, x){
	stopifnot(is.numeric(x), is.numeric(n),
	          length(n) == 1, floor(n) == ceiling(n), n >= 0)
	x <- as.matrix(x)
	N <- nrow(x)
	if (n == 0) return(rep(1, N))
    # generate the Hermite polynomials up to degree n
    Hp <- matrix(0, nrow=nrow(x), ncol=n+1)
    Hp[,1] <- 1; Hp[,2] <- 2*x
    if (n>1){
        for (i in 3:(n+1)) {
	        j <- i-1
	        k <- i-2
	        Hp[,i] <- 2*x*Hp[,j]-2*k*Hp[,k]
        }
    }
    w <- exp(-(x^2)/2)
    c <- 1/sqrt((sqrt(pi)*(2^(n))*factorial(n)))
    hp <- c*Hp[,n+1]*w
    return(hp)
}
#Cosine Polynomial
cospoly <- function(n, x){
	stopifnot(is.numeric(x), is.numeric(n),
	          length(n) == 1, floor(n) == ceiling(n), n >= 0)
	#Normalize x to [0,1]
	x <- (x-min(x))/(max(x)-min(x))
	x <- as.matrix(x)
	return(cos(n*pi*x))
}
#Matrix of Tensor Product Positions
tensormat <- function(D){
	M <- as.numeric(unlist(D))
	#Create Matrix for tensor product coefficient locations
	if (length(unique(M))==1){
		vecs <- mapply(seq, 0, M)
		tmp <- split(vecs, rep(1:ncol(vecs), each=nrow(vecs)))
		tmp <- as.matrix(do.call(expand.grid, tmp))
	} else {
		vecs <- mapply(seq, 0, M)
		tmp <- as.matrix(do.call(expand.grid, vecs))

	}
	return(tmp)
}
tensor <- function(X, D){
	tmp <- tensormat(D=D)
	tensor <- apply(tmp, 1, function(z) cprod(sapply(1:length(z), function(zz) hermite(n=z[zz], x=X[,zz]))))
	return(tensor)
}
tensdens <- function(Xstar, X, C, d, beta){
	return((tensor(D=d, X=cbind(Xstar, X, C))%*%beta)^2)
}
pardens <- function(Y, X, C, d, beta, par){
	#For parametric mean
	if (par==TRUE){
		betamu <- beta[1:(length(beta)-1)]
		betavar <- beta[length(beta)]
		mu <- cbind(1, X, C)
	} else if (par==FALSE){
		# For Semi-parametric mean
		dtot <- prod(unlist(d)+1)
		betamu <- beta[1:dtot]
		betavar <- beta[length(beta)]
		mu <- tensor(D=d, X=cbind(X, C))
	}
	#Log-Normal
	return((-0.5*log(2*pi)-0.5*log(betavar)-0.5/betavar*(Y-(mu%*%betamu))^2))
}
#Score Function for Log-Likelihood
parjac <- function(Y, X, C, d, beta, par){
	#For parametric mean
	if (par==TRUE){
		betamu <- beta[1:(length(beta)-1)]
		betavar <- beta[length(beta)]
		mu <- cbind(1, X, C)
	} else if (par==FALSE){
		#For Semi-parametric mean
		dtot <- prod(unlist(d)+1)
		betamu <- beta[1:dtot]
		betavar <- beta[length(beta)]
		mu <- tensor(D=d, X=cbind(X, C))
	}
	#Jacobian
	jac1 <- -(betavar^-1)*mu*c(Y-(mu%*%betamu))
	jac2 <- -0.5/betavar+(0.5/(betavar^2))*(Y-(mu%*%betamu))^2
	jac <- c(colMeans(jac1), mean(jac2))
	return(jac)
}
hermconst <- function(d, loc, beta){
	if (loc==0){
		c1 <- sum(beta^2)
		return(c1)
	} else if (loc==1){
		dd <- d$dstar
		c1 <- sum(beta^2)
		c2 <- sum(beta[1:dd]*beta[2:(dd+1)]*sqrt(2*(1:dd)))
		return(c(c1, c2))
	}
}
#Function that draws initial values for optimization that satisfies density constraints
hermbinit <- function(d){
	dtot <- prod(unlist(d)+1)
	betainit <- runif(dtot)
	binit <- betainit/sqrt(sum(betainit^2))
	return(binit)
}
fobjz <- function(Xstar, Z, C, dz, beta){
	# dens <- log(tensdens(Xstar=Xstar, X=Z, C=C, d=dz, beta=beta))
	dens <- pardens(Y=Xstar, X=Z, C=C, d=dz, beta=beta, par=TRUE)
	return(-mean(dens))
}
fobjx <- function(Xstar, X, C, dx, beta){
	# dens <- log(tensdens(Xstar=X-Xstar, X=X, C=C, d=dx, beta=beta))
	dens <- pardens(Y=X, X=Xstar, C=C, d=dx, beta=beta, par=TRUE)
	return(-mean(dens))
}
fobjy <- function(Y, Xstar, C, dy, beta){
	densy <- pardens(Y=Y, X=Xstar, C=C, d=dy, beta=beta, par=TRUE)
	return(-mean(densy))
}
dobjy <- function(Y, Xstar, C, dy, beta){
	dobjy <- parjac(Y=Y, X=Xstar, C=C, d=dy, beta=beta, par=TRUE)
	return(dobjy)

}
#Posterior density 
posterior <- function(Y, X, Z, Xstar, C, dy, dx, dz, beta){
	#Position of beta for Y|X,C
	py <- dy
	betay <- beta[1:dy]
	logdensY <- pardens(Y=Y, X=Xstar, C=C, d=dy, beta=betay, par=TRUE)
	#X|Xstar
	# p1x <- py+1; pnx <- py+prod(unlist(dx)+1)
	p1x <- py+1; pnx <- py+dx
	betax <- beta[p1x:pnx]
	# logdensX <- log(tensdens(Xstar=X-Xstar, X=X, C=C, d=dx, beta=betax))
	logdensX <- pardens(Y=X, X=Xstar, C=C, d=dx, beta=betax, par=TRUE)
	#X|Z
	p1z <- pnx+1; pnz <- length(beta)
	betaz <- beta[p1z:pnz]
	# logdensZ <- log(tensdens(Xstar=Xstar, X=Z, C=C, d=dz, beta=betaz))
	logdensZ <- pardens(Y=Xstar, X=Z, C=C, d=dz, beta=betaz, par=TRUE)
	#Combine
	dens <- logdensY+logdensX+logdensZ
	return(dens)

}









