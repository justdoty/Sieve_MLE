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
#Function defining which parameters to constrain
eqfun <- function(X, dz, dx, loc, beta){
	if (is.null(dx)){
		D <- dz
		zloc <- (length(dz$dobs)+1):(length(dz$dobs)+length(dz$dstar))
	} else if (is.null(dz)){
		D <- dx
		zloc <- 1:length(D$dobs)
	}
	M <- tensormat(D=D)
	#Beta0=1
	b0 <- beta[1]
	#Parameter locations for integration constraints
	b1 <- beta[which(rowSums(M[, zloc, drop=FALSE])==0)][-1]
	#Vector for integration contrains
	bint <- c(b0, b1)
	if (loc==1){
		b2loc <- which(cprod(M)==1)
		b21 <- beta[b2loc]
		b22 <- beta[which(cprod(M[, zloc, drop=FALSE])==1)][-b2loc]
		constraints  <- c(bint, b21, b22)
	} else if(loc==0){
		constraints <- bint
	}
	return(constraints)
}
#Function defining the vector for eqfun to be constrained at
eqB <- function(dz, dx, loc){
	if (is.null(dx)){
		D <- dz
		zloc <- (length(dz$dobs)+1):(length(dz$dobs)+length(dz$dstar))
	} else if (is.null(dz)){
		D <- dx
		zloc <- 1:length(D$dobs)
	}
	M <- tensormat(D=D)
	b0 <- 1
	#Parameter locations for integration constraints
	b1intloc <- 1
	b2intloc <- which(rowSums(M[, zloc, drop=FALSE])==0)[-1]
	bintloc <- c(b1intloc, b2intloc)
	#Vector for integration constraint values
	bintval <- c(1, rep(0, length(b2intloc)))
	if (loc==1){
		#For centering restriction(s)
		b1cloc <- which(cprod(M)==1)
		b2cloc <- which(cprod(M[, zloc, drop=FALSE])==1)[-b2loc]
		bcloc <- c(b1cloc, b2cloc)
		#Vector for centering constraint values
		bcval <- c(rep(1, length(b1cloc)), rep(0, length(b2cloc)))
		#Combining the location for parameter constraints
		bloc <- c(bintloc, bcloc)
		#Combining the values for parameter constraints
		bval <- c(bintval, bcval)
	} else if(loc==0){
		#If no normalization (i.e. on f_{xstar|z})
		bloc <- bintloc
		bval <- bintval
	}
	return(list(loc=bloc, val=bval))
}
#Function that draws initial values for optimization that satisfies density constraints
binit <- function(dz, dx, loc){
	if (is.null(dx)){
		D <- dz
		zloc <- (length(dz$dobs)+1):(length(dz$dobs)+length(dz$dstar))
	} else if (is.null(dz)){
		D <- dx
		zloc <- 1:length(D$dobs)
	}
	M <- tensormat(D=D)
	dim <- prod(as.numeric(unlist(D))+1)
	betainit <- runif(dim)
	#Drop values that have constraints
	constr <- eqB(dz=dz, dx=dx, loc=loc)$loc
	binit <- betainit[-constr]
	return(binit)
}
fobj <- function(X, Xstar, dz, dx, loc, beta){
	if (is.null(dx)){
		D <- dz
	} else if (is.null(dz)){
		D <- dx
	}
	#Parameter constraints
	constr <- eqB(dz=dz, dx=dx, loc=loc)
	#Sum over unobservables
	avg <- 0
	for (m in 1:dim(Xstar)[2]){
		xmat <- tensor(D=D, X=cbind(X, Xstar[,m,]))
		#For constrained parameters
		avg1 <- xmat[,constr$loc]%*%constr$val
		#For unconstrained parameters
		avg2 <- xmat[,-constr$loc]%*%beta
		avg <- avg+log((avg1+avg2)^2)
	}
	return(sum(avg))
}
#Jacobian matrix
fjac <- function(X, Xstar, dz, dx, loc, beta){
	if (is.null(dx)){
		D <- dz
	} else if (is.null(dz)){
		D <- dx
	}
	#Parameter constraints
	constr <- eqB(dz=dz, dx=dx, loc=loc)
	#Sum over unobservables
	jac <- 0
	for (m in 1:dim(Xstar)[2]){
		xmat <- tensor(D=D, X=cbind(X, Xstar[,m,]))
		#For constrained parameters
		x1 <- xmat[,constr$loc]%*%constr$val
		#For unconstrained parameters
		x2 <- xmat[,-constr$loc]%*%beta
		xtot <- x1+x2
		jac <- jac+colSums(2*xmat[,-constr$loc]/c(xtot))
		
	}
	return(jac)

}
#Posterior density 
posterior <- function(X, Xstar, dz, dx, beta, loc){
	if (is.null(dx)){
		D <- dz
	} else if (is.null(dz)){
		D <- dx
	}
	xmat <- tensor(D=D, X=cbind(X, Xstar))
	#Parameter constraints
	constr <- eqB(dz=dz, dx=dx, loc=loc)
	#For constrained parameters
	x1 <- xmat[,constr$loc]%*%constr$val
	#For unconstrained parameters
	x2 <- xmat[,-constr$loc]%*%beta
	dens <- log((x1+x2)^2)
	dens <- data.table(idvar=idvar, dens=dens)[,sum(dens),keyby=idvar]$V1
	return(dens)

}









