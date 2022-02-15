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
    c <- sqrt(1/sqrt(pi)*factorial(n+1)*2^(n+1))
    hp <- c*Hp[,n+1]*w
    return(hp)
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
tensor <- function(X, D, basis){
	tmp <- tensormat(D=D)
	tensor <- apply(tmp, 1, function(z) cprod(sapply(1:length(z), function(zz) legendre(n=z[zz], x=X[,zz]))))
	return(tensor)
}
# #Function defining the vector for eqfun to be constrained at
# eqB <- function(dz, dx, loc){
# 	if (is.null(dx)){
# 		D <- dz
# 		zloc <- (length(D$dobs)+1):(length(D$dobs)+length(D$dstar))
# 	} else if (is.null(dz)){
# 		D <- dx
# 		zloc <- 1:length(D$dobs)
# 	}
# 	M <- tensormat(D=D)
# 	b0 <- 1
# 	#Parameter locations for integration constraints
# 	b1intloc <- 1
# 	b2intloc <- which(rowSums(M[, zloc, drop=FALSE])==0)[-1]
# 	bintloc <- c(b1intloc, b2intloc)
# 	#Vector for integration constraint values
# 	bintval <- c(1, rep(0, length(b2intloc)))
# 	if (loc==1){
# 		#For centering restriction(s)
# 		b1cloc <- which(cprod(M)==1)
# 		b2cloc <- which(cprod(M[, zloc, drop=FALSE])==1)[-b2loc]
# 		bcloc <- c(b1cloc, b2cloc)
# 		#Vector for centering constraint values
# 		bcval <- c(rep(1, length(b1cloc)), rep(0, length(b2cloc)))
# 		#Combining the location for parameter constraints
# 		bloc <- c(bintloc, bcloc)
# 		#Combining the values for parameter constraints
# 		bval <- c(bintval, bcval)
# 	} else if(loc==0){
# 		#If no normalization (i.e. on f_{xstar|z})
# 		bloc <- bintloc
# 		bval <- bintval
# 	}
# 	bconstr <- list(loc=bloc, val=bval)
# 	return(bconstr)
# }
#Function that draws initial values for optimization that satisfies density constraints
binit <- function(dy, dz, dx){
	if (is.list(dy)){
		dtot <- prod(unlist(dy)+1)+prod(unlist(dx)+1)+prod(unlist(dz)+1)
	} else {
		dtot <- dy+prod(unlist(dx)+1)+prod(unlist(dz)+1)
	}
	betainit <- runif(dtot)
	binit <- betainit
	return(binit)
}
fobjz <- function(Xstar, Z, C, dz, beta){
	densz <- 0
	for (m in 1:dim(Xstar)[2]){
		densz <- densz+log((tensor(D=dz, X=cbind(Z, Xstar[,m,], C))%*%beta)^2)
	}
	densz <- densz/dim(Xstar)[2]
	return(mean(densz))
}
dobjz <- function(Xstar, Z, C, dz, beta){
	densz <- 0
	for (m in 1:dim(Xstar)[2]){
		tensz <- tensor(D=dz, X=cbind(Z, Xstar[,m,], C))
		densz <- densz+colMeans(2*tensz/c(tensz%*%beta))
	}
	densz <- densz/dim(Xstar)[2]
	return(densz)

}
fobjx <- function(Xstar, X, C, dx, beta){
	densx <- 0
	for (m in 1:dim(Xstar)[2]){
		densx <- densx+log((tensor(D=dx, X=cbind(X, Xstar[,m,], C))%*%beta)^2)
	}
	densx <- densx/dim(Xstar)[2]
	return(mean(densx))
}
dobjx <- function(Xstar, X, C, dx, beta){
	densx <- 0
	for (m in 1:dim(Xstar)[2]){
		tensx <- tensor(D=dx, X=cbind(X, Xstar[,m,], C))
		densx <- densx+colMeans(2*tensx/c(tensx%*%beta))
	}
	densx <- densx/dim(Xstar)[2]
	return(densx)

}
fobjy <- function(Y, Xstar, C, dy, densY, beta){
	densy <- 0
	for (m in 1:dim(Xstar)[2]){
		if (missing(densY)){
			densy <- densy+log((tensor(D=dy, X=cbind(Y, Xstar[,m,], C))%*%beta)^2)
		} else{
			densy <- densy+densY(Y=Y, Xstar=Xstar[,m,], C=C, beta=beta)
		}	
	}
	densy <- densy/dim(Xstar)[2]
	return(-mean(densy))
}
dobjy <- function(Y, Xstar, C, dy, jacY, beta){
	dobjy <- 0
	for (m in 1:dim(Xstar)[2]){
		if (missing(jacY)){
			tensy <- tensor(D=dy, X=cbind(Y, Xstar[,m,], C))
			dobjy <- dobjy+colMeans(2*tensy/c(tensy%*%beta))
		} else {
			dobjy <- dobjy+jacY(Y=Y, Xstar=Xstar[,m,], C=C, beta=beta)
		}
	}
	dobjy <- dobjy/dim(Xstar)[2]
	return(dobjy)

}
#Posterior density 
posterior <- function(densY, Y, X, Z, Xstar, C, dy, dx, dz, beta){
	#Position of beta for Y|X,C
	if (missing(densY)){
		py <- prod(unlist(dy)+1)
		betay <- beta[1:py]
		densY <- tensor(D=dy, X=cbind(Y, Xstar, C))%*%betay
		densY <- log(densY^2)
	} else {
		py <- dy
		betay <- beta[1:dy]
		densY <- densY(Y=Y, Xstar=Xstar, C=C, beta=betay)
	}
	#Position of beta for X|Xstar, C
	p1x <- py+1; pnx <- py+prod(unlist(dx)+1)
	betax <- beta[p1x:pnx]
	#Position of beta for Xstar|Z, C
	p1z <- pnx+1; pnz <- length(beta)
	betaz <- beta[p1z:pnz]
	#Densities
	densx <- tensor(D=dx, X=cbind(X, Xstar, C))%*%betax
	densz <- tensor(D=dz, X=cbind(Z, Xstar, C))%*%betaz
	#Combine
	dens <- densY+log(densx^2)+log(densz^2)
	dens <- data.table(id=idvar, dens=as.numeric(dens))[,sum(dens),keyby=id]$V1
	return(dens)

}









