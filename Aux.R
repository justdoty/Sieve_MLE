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
	tensor <- apply(tmp, 1, function(z) cprod(sapply(1:length(z), function(zz) legendre(n=z[zz], x=X[,zz]))))
	return(tensor)
}
#Constraint Vector
cvec <- function(dz, dx, loc, beta){
	if (missing(dx)){
		D <- dz
		zloc <- (length(dz$dobs)+1):(length(dz$dobs)+length(dz$dstar))
	} else if (missing(dz)){
		D <- dx
		zloc <- 1:length(D$dobs)
	}
	M <- tensormat(D=D)
	#Beta0=1
	b0 <- beta[1]-1
	#Parameter locations for integration constraints
	b1 <- beta[which(rowSums(M[, zloc, drop=FALSE])==0)][-1]
	#Vector for integration contrains
	bint <- c(b0, b1)
	if (loc==1){
		b2loc <- which(cprod(M)==1)
		b21 <- beta[b2loc]-1
		b22 <- beta[which(cprod(M[, zloc, drop=FALSE])==1)][-b2loc]
		bloc <- c(bint, b21, b22)
	} else if(loc==0){
		bloc <- bint
	}
	return(bloc)
}
#Function that draws initial values for optimization that satisfies density constraints
binit <- function(dz, dx, loc){
	if (missing(dx)){
		D <- dz
		zloc <- (length(dz$dobs)+1):(length(dz$dobs)+length(dz$dstar))
	} else if (missing(dz)){
		D <- dx
		zloc <- 1:length(D$dobs)
	}
	M <- tensormat(D=D)
	dim <- prod(as.numeric(unlist(D))+1)
	betainit <- runif(dim)
	#Beta0=1
	betainit[1] <- 1
	#Parameter locations for integration constraints
	betainit[which(rowSums(M[, zloc, drop=FALSE])==0)][-1] <- 0
	if (loc==1){
		b2loc <- which(cprod(M)==1)
		betainit[b2loc] <- 1
		betainit[which(cprod(M[, zloc, drop=FALSE])==1)][-b2loc] <- 0
		binit <- betainit
	} else if(loc==0){
		binit <- betainit
	}
	return(binit)
}
fobj <- function(X, D, dimun, beta){
	#Draw values of unobservables from posterior distribution
	Xstar <- MCMC(X=X, D=D, dimun=dimun, bdraws=bdraws, mdraws=mdraws, beta=beta)
	#Sum over unobservables
	avg <- 0
	for (m in 1:dim(Xstar)[2]){
		avg <- avg+log(tensor(J=J, X=cbind(X, Xstar[,m,]))%*%beta)
	}
	return(sum(avg))
}
#Posterior density of Xstar|Z
posterior <- function(X, Xstar, D, beta){
	dens <- log(tensor(D=D, X=cbind(X, Xstar))%*%beta)
	dens <- data.table(id=idvar, dens=dens)[,sum(dens),keyby=id]$V1
	return(dens)

}









