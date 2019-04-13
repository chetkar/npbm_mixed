
 fn.offset <- function(parm){

	for(xx  in 1:parm$p ){

		theta.v      <- parm$clust$A.mt[,parm$clust$c.v[xx]]
		post.prec <- sum( (1/parm$clust$sdm.mt[,xx])^2) + 1
	    post.sd   <- sqrt(1/post.prec)
		
		post.mean <- (post.sd)^2*sum( (parm$X[,xx])/parm$clust$sdm.mt[,xx]^2)
		
			
		prop1     <- rnorm(1, post.mean, post.sd)
				
		if( parm$data.type[xx] == 1){
			parm$offset[xx] <- mean( qnorm( ( parm$OX[,xx] + 0.01 )/1.02 ) )
		}
		
		if( parm$data.type[xx]  == 4 ){
			parm$offset[xx] <- mean( qnorm( ( parm$OX[,xx] + 0.01 )/5.02 ) ) 
		}
		
		if( parm$data.type[xx]  == 3 ){
			 parm$OX   <- as.matrix(parm$OX)
			 parm$offset[xx] <- log(mean(parm$OX[,xx]) + 0.01)
		}
		
		if( parm$data.type[xx] == 5 ){
			 parm$offset[xx] <- prop1
		}
		
		if( parm$data.type[xx] == 2  ){
			parm$offset[xx] <- prop1
		}
		
	 }
		 
	 for( xx in 1:parm$p ){
	 		parm$X[,xx] <- parm$X[,xx] - parm$offset[xx]
		}
		
 parm	
 }

fn.transform <- function(parm, data){

  # data$OX
  n  <- dim(data$OX)[1]
  p  <- dim(data$OX)[2]
  TOX <- TX <- matrix(NA,n,p)
  
  
 for( xx in 1:p) {	

	theta.v      <- parm$clust$A.mt[,parm$clust$c.v[xx]] + parm$offset[xx]
 
	#Check logit	
	if(parm$data.type[xx] == 1)
	{
	 ind0      <- which(data$OX[,xx] == 0)	
	 ind1      <- which(data$OX[,xx] == 1)
	
	 TX[ind0,xx] 	  <- rtnorm( length(ind0) ,mean = theta.v[ind0] , sd = 1, upper = 0)
	 TX[ind1,xx] 	  <- rtnorm( length(ind1) ,mean = theta.v[ind1] , sd = 1, lower = 0)	
	
	} 
	
	# Check Poisson
	if(parm$data.type[xx] == 3)
	{
	TX[,xx] 	  	  <- theta.v + exp(-theta.v)*(data$OX[,xx] ) - 1
	} 
		
	#Check Ordinal
	if(parm$data.type[xx] == 4)
	{
    
	v  <- sort(unique(data$OX[,xx]))
	new.cutoff <- cutoff  <- parm$data.cutoff[,xx]
	h   <- min(data$OX[,xx])
			
	if(h != 5 ){
			hh <- h + 1
			#One can fix this to be -1?
			new.cutoff[hh] <- cutoff[hh] <- -1 
			new.cutoff[h]  <- cutoff[h]  <- -Inf
	}
	
	h   <- max(data$OX[,xx])
	hh  <- h + 1
	new.cutoff[hh]  <- cutoff[hh]  <- +Inf
	    	
	for(j in 1: length(v) )
	{
		jj  <- j + 1
		ind <- which(data$OX[,xx] == v[j])
		
		
		if( j != length(v) ){
			TX[ind,xx] <- rtnorm(length(ind), mean = theta.v[ind], sd = 1, lower= cutoff[v[j]], upper= cutoff[v[jj]] )
		}else{
			TX[ind,xx] <- rtnorm(length(ind), mean = theta.v[ind], sd = 1, lower= cutoff[v[j]], upper= cutoff[v[length(v)] + 1])
		}
	}
	
	
	ll <- length(v) -1
	
	for(j in 2: length(v)-1){
		ind   <- which(data$OX[,xx] == v[j] )
		jj    <- j + 1
		 
			ind.u <- which(data$OX[,xx] == v[jj] )
			
			if(length(ind) > 0 ){
				max.i <- max(TX[ind,xx])
			}else{
				max.i <- cutoff[v[j]]
			}
    	
		 if(length(ind.u) > 0 ){
		     min.i <- min(TX[ind.u,xx])
		  }else{
			 jj <- j + 2
			
			if(jj == length(v) + 1 ){
				min.i <- cutoff[v[length(v)] + 1] 
			}else{
				min.i <- cutoff[v[jj]]
			} 
		 
		 }
	
		  lower <- max( max.i,cutoff[v[j]]) 
		  jj <- j + 2
		  
		  if(jj == length(v) + 1 ){
			upper <- min( min.i,cutoff[v[length(v)] + 1]) 
		  }else{
			upper <- min( min.i,cutoff[v[jj]])
		  }
		  
		  jj <- j + 1
		
		new.cutoff[v[jj]] <- runif(1,lower, upper)
		}
	
		

	parm$data.cutoff[,xx] <- new.cutoff 
	}
	
	if(parm$data.type[xx] == 5 ){
		
		#Numerical Stability			
		z <- data$OX[,xx]
		
		ind0 <- which(z < 0.005 )
		ind1 <- which(z > 0.995 )
		z[ind0] = z[ind0] + 0.005
		z[ind1] = z[ind1] - 0.005
		
		y    <- log( z/(1- z) )
		TX[,xx] <- y		

	}
	
	if(parm$data.type[xx] == 2 ){
 
		TX[,xx] <- data$OX[,xx]
	}
		
}


parm$X  <- TX

parm	
}


#Called up by elementwise_DP.functions.R
fn.matrix.sd <- function(parm,data){

 n <- nrow(parm$X)
 p <- ncol(parm$X)
 parm$clust$sd.matrix.mt <- NULL
 parm$clust$sd.matrix.mt <- matrix(0,n,p)
 parm$clust$s.full.mt   <- matrix(0,n,p)


 
 for(l in 1:p ){

	theta.v <- parm$clust$A.mt[,parm$clust$c.v[l] ]
	s.v     <- parm$clust$s.mt[,parm$clust$c.v[l]]
	parm$clust$s.full.mt[,l] <- s.v
  
 }# end forloop

 
 parm 

}


# Called by PDP_fn.log.lik
fn.equivsd <- function(parm,data){

 n <- nrow(parm$X)
 p <- ncol(parm$X)
 parm$clust$sdm.mt <- NULL
 parm$clust$sdm.mt <- matrix(NA,n,p)


  
 for(l in 1:p ){

	theta.v <- parm$clust$A.mt[,parm$clust$c.v[l]] + parm$offset[l]

	if ( parm$data.type[l] == 1){
		
		sd.v.m  <- rep(1, length(theta.v) )

	  }
	  

	if ( parm$data.type[l] == 2 | parm$data.type[l] == 6){
		
		 sd.v.m  <- rep(parm$tau, length(theta.v) )
	  }
	  
	if ( parm$data.type[l] == 3){
		
 		 sd.v.m  <- sqrt( exp(-theta.v) )
	  }  

	if ( parm$data.type[l] == 4){
		
 		 sd.v.m  <- rep(1,length(theta.v) )
	  }
	
	if (parm$data.type[l] == 5){
		
			
		a1   <- exp(theta.v)
		a2   <- ( 1 + exp(theta.v) )
		mu   <- a1/a2
		
		
		#Numerical Stability
		z    <- mu
		
		phi  <- parm$phi
		mu_new <- digamma(mu*phi)  - digamma( (1-mu)*phi )
				
		phi  <- parm$phi
		mu_new <- digamma(mu*phi)  - digamma( (1-mu)*phi )
				
		var1 <- trigamma(mu*phi) + trigamma( (1-mu)*phi)
		var2 <- (mu*(1-mu))^2
		
		sd.v.m <- sqrt(var1 * var2)	
			  
		}
	

  parm$clust$sdm.mt[,l] <- sd.v.m

 }# end forloop

 
 parm$clust$max.sd <- max(parm$clust$sdm.mt)
 parm$clust$min.sd <- min(parm$clust$sdm.mt)
 parm$clust$avg.sd <- mean(parm$clust$sdm.mt)
 
 parm 

}


# Called by PDP_fn.log.lik
fn.margsd <- function(parm,I.k){

 for(l in 1:length(I.k) ){

	theta.v <- parm$clust$A.mt[,parm$clust$c.v[I.k[l]]] + parm$offset[I.k[l]]

	if ( parm$data.type[l] == 1){
		
		sd.v.m  <- rep(1, length(theta.v) )
	  }
	  

	if ( parm$data.type[l] == 2 | parm$data.type[l] == 6){
		
		 sd.v.m  <- rep(parm$tau, length(theta.v) )
	  }
	  
	if ( parm$data.type[l] == 3){
		
 		 sd.v.m  <- sqrt( exp(-theta.v) )
	  }  

	if ( parm$data.type[l] == 4){
		
 		 sd.v.m  <- rep(1,length(theta.v) )
	  }
	
	if (parm$data.type[l] == 5){
		
			
		a1   <- exp(theta.v)
		a2   <- ( 1 + exp(theta.v) )
		mu   <- a1/a2
				
		#Numerical Stability
		z    <- mu
		
		phi  <- parm$phi
		mu_new <- digamma(mu*phi)  - digamma( (1-mu)*phi )
				
		var1 <- trigamma(mu*phi) + trigamma( (1-mu)*phi)
		var2 <- (mu*(1-mu))^2
		
		sd.v.m <- sqrt(var1 * var2)	
			  
		}
	
  parm$clust$sdm.mt[,I.k[l]] <- sd.v.m

 }# end forloop

 
 parm$clust$max.sd <- max(parm$clust$sdm.mt)
 parm$clust$min.sd <- min(parm$clust$sdm.mt)
 parm$clust$avg.sd <- mean(parm$clust$sdm.mt)
 
 parm 

}
post.phi <- function(phi,parm,data)
    {	
	   log.fn  <- 0
	   rate.parm <- 0.01
	   shape.parm <- 0.01
	
	   for(l in 1:parm$p){		
		 if(parm$data.type[l] == 5){
			   theta.v <- parm$clust$A.mt[,parm$clust$c.v[l]]	
			   a1   <- exp(theta.v)
			   a2   <- ( 1 + exp(theta.v) )
		       mu   <- a1/a2
			   			   
				#Numerical Stability
				z    <- mu
				
				var1 <- trigamma(mu*phi) + trigamma( (1-mu)*phi)
				var2 <- (mu*(1-mu))^2
				sd   <- sqrt(var1*var2)
				log.fn <- log.fn + sum( dnorm(parm$X[,l], theta.v, sd, log =TRUE))
		  }
	   }
	   log.prior <- dgamma(phi, shape.parm, rate.parm , log = TRUE)
	   log.post <- log.prior + log.fn
	
    log.post
   }

fn.sample.phi <- function(parm, data){
			
			phi0 	 <- rgamma(1, parm$phi , 1) + 0.01
			phi1 	 <- parm$phi
			log.rat <- dgamma(phi1, phi0, 1, log = TRUE) - dgamma(phi0, phi1 , 1, log =TRUE) 	
			log.rat <- log.rat + max(min(post.phi(phi0,parm,data),10^10),-10^10) - max(min(post.phi(phi1,parm,data),10^10),-10^10) 
		
			rat     <- exp(log.rat)
			uni     <- runif(1,0,1)
	
		 if( uni < rat ){
			 parm$phi <- phi0
			}

    parm	
}

fn.sample.M <- function(parm){
#Samples from a posterior density of Mixture Gamma Distribution
# Paper Reference Escobar and West (1995)
# a, b are hyperparameters
    a <- 0.001
	b <- 0.001
     
	first.parm <- parm$clust$M + 1
	sec.parm   <- parm$n2*parm$clust$G

	nita <- rbeta(1, first.parm, sec.parm )

	p   <- (a + parm$clust$K - 1)/( parm$n2*(b - log(nita) ) + a + parm$clust$K - 1)
	ind <- rbinom(1, 1, p)
	
	if(ind == 1){
	a  <- rgamma(1, shape = a + parm$clust$K, rate= b - log(nita) )
	}else{
	a  <- rgamma(1, shape = a + parm$clust$K - 1, rate = b - log(nita) )	
	}

	parm$clust$M <- a

parm	
}



fn.dmvnorm <- function(x, mean, sigma, inv.sigma, log=TRUE)
	{

	# Computes multivariate normal density function
	# a little faster than dmvnorm function of R!

	if (missing(inv.sigma))
		{inv.sigma <- solve(sigma)
		}

	logdet <- as.numeric(determinant(inv.sigma, logarithm=TRUE)$mod)
	r <- length(x)
	Q <- colSums(inv.sigma * (x-mean))
	Q <- sum(Q * (x-mean))

	val <- -r/2*log(2*pi) + logdet/2 - Q/2
	if (!log)
		{val <- exp(val)
		}

	val
	}



fn.quality.check <- function(parm)
	{err <- 0

	if (parm$tBB_flag)
	  {
	  if (sum(diag(parm$clust$tBB.mt) < 0) > 0)
		  {err <- 2
	    }
	  }

	if (ncol(parm$clust$A.mt) != parm$clust$G)
		{err <- 3
		}

	if (ncol(parm$clust$B.mt) != (parm$clust$G+1))
		{err <- 4
		}

	if ((sum(parm$clust$C.m.vec) + parm$clust$C.m0) != parm$p)
		{err <- 5
		}

	if (length(parm$clust$n.vec) != parm$clust$K)
		{err <- 6
		}

	if (length(parm$clust$phi.v) != parm$clust$K)
		{err <- 7
		}

	if ((sum(parm$clust$n.vec) + parm$clust$n0) != parm$N)
		{err <- 8
		}

	err
	}

######################

fn.init.clusters <- function(parm)
{
		
	# Not Working for the discrete case
    n <- dim(parm$X)[1]
	p <- dim(parm$X)[2]

	X.mt  <- matrix(as.vector(parm$X),n,p)
	
	
	num.centers <- parm$G.new		

	options(warn=2)
	
	tmp2 <- kmeans( t( X.mt ), iter.max=1000, centers=num.centers, nstart=2)
	
	parm$clust$c.v <- tmp2$cluster
	parm$clust$G <- length(tmp2$size)
	

	parm$clust$C.m.vec <- array(,parm$clust$G)	
		
		for (g in 1:parm$clust$G)
		{I.g <- (parm$clust$c.v==g)
		 parm$clust$C.m.vec[g] <- sum(I.g)
		}	
	
	

	#
	
	###########################

	# start from PDP model with d=0 (i.e DP)
	#parm$d <- 1/3
	parm$d <- 0

	parm
}


fn.eda <- function(parm, data, computeMode)
{

	parm <- fn.init.clusters(parm)
	# reintroduced on 6/29/12
	parm$G.max <- min(parm$p/2, round(parm$clust$G*1.1))
	parm <- fn.poissonDP.hyperparm(data, parm, w=.01, max.d=1)

	parm$Y <- parm$clust$A.mt <- array(,c(parm$n2,parm$clust$G))
	parm$clust$C.m.vec <- array(,parm$clust$G)

	#Probit or Ordinal (on least priority as it isn't scale dependent)
	if(sum(parm$data.type== 1) > 0 ){
	   parm$clust$mu2   <- 0
	   ind <- which(parm$data.type == 1)	
	   parm$clust$tau2 <- 1
	}
	
	#Divided by number of categories
	if(sum(parm$data.type== 4) > 0 ){
	   parm$clust$mu2  <- 0
	   parm$clust$tau2 <- 1
	}
	#Normal
	if(sum(parm$data.type == 2) > 0 ){
	   ind <- which(parm$data.type == 2)
	   parm$offset[ind]  <- mean(as.vector(parm$X[,ind]))
	   parm$clust$mu2  <- 0
	   parm$clust$tau2 <- diff(range(as.vector(parm$X[,ind])))/6
	}
	
	#Poisson
	if(sum(parm$data.type == 3) > 0 ){
	  ind <- which(parm$data.type == 3)	
	   parm$clust$mu2  <- 0
	   parm$offset[ind] <- 0
	   parm$clust$tau2 <-   log( diff(range(as.vector(parm$X[,ind])))/3)	
	   
	}

	#Proportion
	if( sum(parm$data.type == 5) > 0 ){
	  ind <- which(parm$data.type == 5)	
	  parm$clust$mu2  <-  log( max( mean(as.vector(parm$X[,ind]))/(1 - mean(as.vector(parm$X[,ind]))) , 0.5)  )
	  parm$clust$tau2 <-  diff(range(as.vector(parm$X[,ind])))/2  
	}		
	
		for (g in 1:parm$clust$G)
		{I.g <- (parm$clust$c.v==g)
		 parm$clust$C.m.vec[g] <- m.g <- sum(I.g)
		x.g.v <- parm$X[,I.g]
		 if (m.g > 1)
			{x.g.v <- rowMeans(x.g.v)
			}
	
			parm$Y[,g] <- x.g.v
			
		}

	parm$clust$C.m0 <- parm$p - sum(parm$clust$C.m.vec)

	parm$clust$M <- parm$a.R
	parm$clust$M0 <- .01*parm$clust$M
	
	parm$clust$K <- data$K.max	

    
		
	
	#################################

	parm$g <- rep(1:parm$clust$G,each=parm$n2)
	parm$N <- parm$clust$G * parm$n2

	parm$Y <- as.vector(parm$Y)

	# if (computeMode$useR) {

	tmp <- tryCatch({
	  iter.max <- ifelse((length(parm$Y) > 1000), 10, 1000)
	  kmeans(parm$Y, iter.max = iter.max, centers = data$K.max, nstart = 10,
	         algorithm = "Hartigan-Wong" # TODO: MAcQueen works better?
	  )}, error = function(e) {
	    print("Kmeans did not converge ... using random assignment")
	    cluster <- sample(1:data$K.max, size = length(parm$Y), replace = TRUE)
	    centers <- sapply(1:data$K.max, FUN = function(x) {
	      mean(parm$Y[which(cluster == x)])
	    })
	    list(cluster = cluster, centers = centers, size = data$K.max)
	  })


	parm$clust$s.v <- tmp$cluster
	parm$clust$phi.v <- as.vector(tmp$centers)

	if( sum( parm$data.type == 3 ) > 0 ){
	
		 ind0 <- which(parm$clust$phi.v > 0 )
		 ind1 <- which(parm$clust$phi.v < 0 )
		parm$clust$phi.v[ind0] <- log(parm$clust$phi.v[ind0])
		parm$clust$phi.v[ind1] <- log(-parm$clust$phi.v[ind1])
	
	}	

	parm$clust$n.vec <- tmp$size
	# number of s equal to 0
	parm$clust$n0 <- 0

	parm$clust$s.mt <- array(parm$clust$s.v, c(parm$n2,parm$clust$G))

	for (g in 1:parm$clust$G)
		{parm$clust$A.mt[,g] <- parm$clust$phi.v[parm$clust$s.mt[,g]]
		}

	sum.resid.sq <- 0
	
	
	for (g in 1:parm$clust$G)
		{flag.v <- parm$clust$c.v == g & parm$data.type ==2
		X.g.mt <- parm$X[,flag.v]
		a.g.v <- parm$clust$A.mt[,g]
		resid.g.mt <- X.g.mt - a.g.v
		sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
		}
	
	p1 <- sum(parm$data.type == 2)
	parm$tau_int <- parm$tau <- sqrt(sum.resid.sq/parm$n2/p1)

	###################################

	parm$tau_0 <- sqrt(1+parm$tau^2)


	## objects of full size (based on all n2 cases)
	parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)

	if (parm$tBB_flag)
	  {parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt
	  }

	parm <- fn.assign.priors(parm, data)

  parm

	}



fn.gen.clust <- function(parm, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, computeMode)
	{

  ###########################################
  # Missing X values
  ###########################################

   n <- dim(data$X)[1]
   p <- dim(data$X)[2]

   X.mt <- data$X
	
  parm$X <- X.mt
  parm$num.X.miss <- sum(is.na(parm$X))
  tmp <- which(is.na(parm$X), arr=TRUE)
  parm$X.missing.x <- tmp[,1]
  parm$X.missing.y <- tmp[,2]

 parm$X <- matrix(as.numeric(unlist(parm$X)),nrow(parm$X),ncol(parm$X))  
	# Impute any missing X values by their column-specific means
	# + a small error term to guarantee non-tied values

 	tmp.mean.v <- apply(parm$X, 2, median, na.rm=TRUE)
	tmp.sd.v <- apply(parm$X, 2, sd, na.rm=TRUE)
	if (parm$num.X.miss>0)
	  { 	for (j in 1:parm$p)
		      {indx.j <- is.na(parm$X[,j])
		      if (sum(indx.j) > 0)
			      {parm$X[indx.j,j] <- tmp.mean.v[j] + rnorm(n=sum(indx.j), sd=tmp.sd.v[j]/5)
		      }
	  }
	}

	##################

	parm$G.new <- data$G.max
	#print(2)
	parm$offset <- rep(0,parm$p)
	parm <- fn.eda(parm, data, computeMode)

	#################
	#print(3)
	parm <- fn.hyperparameters(data, parm)
	#parm <- fn.sample.bet(parm,data)
	#Initializing parm$phi
		
	parm$phi    <- 10


	parm <- fn.equivsd(parm,data)
	parm <- fn.matrix.sd(parm,data)


	 if( sum(parm$data.type == 3) > 0 ){
	
	 ind0 <- which(parm$X > 0 )
	 ind1 <- which(parm$X < 0 )
	 parm$X[ind0] <- log(parm$X[ind0])
	 parm$X[ind1] <- log(-parm$X[ind1])
	}
	
	parm <- fn.element.DP(data, parm, max.row.nbhd.size, row.frac.probes, computeMode)
	#print(6)
	parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)

	if (parm$tBB_flag)
	  {parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt
	  }

	parm

	}

fn.init <- function(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm, tBB_flag, standardize.X, flip.sign, computeMode = "R" )
	{
	#print(1)

#	parm <- true_parm

	parm <- NULL
	parm$OX <- data$OX

	parm$data.type <- data$data.type
		
	if(sum(parm$data.type == 4) > 0 )
	{
		ind <- which(parm$data.type == 4)
		
		l   <- max(parm$OX[,ind])
		ll  <- l + 1
		ext.cut <- seq(0,l-2,1) - 1
		ext.cut <- c( -Inf, ext.cut, Inf)
		ext.cut[2] <- -1	
		
		parm$data.cutoff <- matrix( ext.cut, ll , ncol(data$OX))
		
	}
    	

	parm$tBB_flag <- tBB_flag
	parm$standardize.X <- standardize.X
	parm$flip.sign <- flip.sign

	parm$n2 <- dim(data$X)[1] # TODO Check
	parm$p <- dim(data$X)[2]  # TODO Check

	### ASSUMING POSITIVE ORIENTATION FOR ALL PDP CLUSTERS
	### IN INITIALIZATION
	parm$clust$orient.v <- rep(1,parm$p)

	# mass parameter of elementwise(s) groups
	# stored later in parm$clust$M
	parm$a.R <- true$a.R

	# mass paramater of columns
	parm$b1 <- true$b1

	# mass paramater of column-intercept cluster
	parm$b0 <- true$b0


	############################
	# For delta neighborhoods
	############################

	parm$col.delta <- .05

	# delta-neighborhood threshold for elements
	parm$row.delta <- .1

	#########################################
	# generating the R- and C- clusters
	########################################

	parm$shift <- true$shift
	
	parm <- fn.gen.clust(parm, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, computeMode)

	parm <- fn.assign.priors(parm, data)

	parm

	}


fn.gen.missing.X <- function(data, parm)
{
  # impute missing X values
  #  X.mt <- data$X
  # Discrete Case		
    #Data <- fn.transform(parm,data)
    X.mt <- data$X
	
  if (parm$num.X.miss > 0)
  {
    for (cc in 1:parm$num.X.miss)
    {i.cc <- parm$X.missing.x[cc]
    j.cc <- parm$X.missing.y[cc]
    c.cc <- parm$clust$c.v[j.cc]
    if (c.cc != 0)
    {mean.cc <- parm$clust$A.mt[i.cc, c.cc]
    }
    if (c.cc == 0)
    {mean.cc <- 1
    }
    X.mt[i.cc, j.cc] <- rnorm(n=1, mean=mean.cc, sd=parm$tau)
    }
  }
  parm$X <- X.mt

  parm
}


fn.standardize_orient.X <- function(parm)
{

  ####
  ## STANDARDIZE X columns to unit variance and zero mean
  #####
  # Do only for columns with NA's
  # For other columns, it's just a one-time calculation at the beginning of MCMC

  if (parm$num.X.miss > 0)
    {tmp.X <- matrix(parm$X[,parm$X.missing.y],col=parm$num.X.miss)
    mean.v <- colMeans(tmp.X)
    sd.v <- apply(tmp.X, 2, sd)
    parm$X[,parm$X.missing.y] <- t((t(tmp.X) - mean.v)/sd.v)
   }

  ####
  ## ORIENT X
  ####

  parm$X <- t(t(parm$X) * parm$clust$orient.v)

  parm
}


fn.assign.priors <- function(parm, data)
	{

	parm$prior$tau <- NULL
	parm$prior$tau$alpha.tau <- 1e-2
	parm$prior$tau$beta.tau <- 1e-2

	parm$prior$tau$max <- sqrt(.75)*sd(as.vector(as.numeric(as.matrix(data$X))), na.rm=TRUE)
	parm$prior$tau$min <- 1e-10
	parm$prior$tau.sq$max <- parm$prior$tau$max^2
	parm$prior$tau.sq$min <- parm$prior$tau$min^2
	parm$prior$inv.tau.sq$max <- 1/parm$prior$tau.sq$min
	parm$prior$inv.tau.sq$min <- 1/parm$prior$tau.sq$max
	
	parm$prior$mu2$mean <- 0
	parm$prior$mu2$sd   <- 1
	parm$prior$tau2$alpha <- 0.01
	parm$prior$tau2$beta <- 0.01
	
	parm$prior$mu$mean <- 0
	parm$prior$mu$sd   <- 1

	parm
	}

#######################################
###
########################################

fn.sample.mu <- function(parm,data){

 ind <- which(parm$data.type == 2)
 x   <- parm$OX[,ind]
 a   <- parm$clust$A.mt[,parm$clust$c.v[ind]]
 
 post.prec <- (1/parm$prior$mu$sd^2 + length(x)/parm$tau^2)
 post.mean <- (1/post.prec)*( parm$prior$mu$mean/parm$prior$mu$sd^2 + sum(x-a)/parm$tau^2 ) 
 post.sd   <- sqrt(1/post.prec)
 
 parm$mu <- rnorm(1, post.mean, sd = post.sd)
 
 parm
 }
	
########################################
#### This is for the continuous case ###
########################################

fn.gen.tau  <- function(data, parm)
	{
	###################
	# update tau
	###################

	# only covariates assigned to non-zero row and non-zero group clusters matter
	# Need to update the clusters to only continuous columns
  #if(!parm$discrete){

	sum.resid.sq <- 0
	count <- 0


	for (g in 1:parm$clust$G)
		{flag.v <- ( (parm$clust$c.v == g & parm$data.type == 2) | (parm$clust$c.v == g & parm$data.type == 6 ) ) 
		z.g.v   <- parm$clust$s.mt[,g] > 0
		
		if ((sum(z.g.v) > 0) & (sum(flag.v)>0) )
			{
		
			X.g.mt <- parm$X[z.g.v,flag.v]
			a.g.v <- parm$clust$A.mt[z.g.v,g]
			resid.g.mt <- X.g.mt - a.g.v
			sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
			count <- count + sum(z.g.v)*sum(flag.v)
			}

		}
	
	shape <- parm$prior$tau$alpha + count/2
	rate <- parm$prior$tau$beta + sum.resid.sq/2

	u.min <- pgamma(parm$prior$inv.tau.sq$min,shape=shape, rate=rate)
	u.max <- pgamma(parm$prior$inv.tau.sq$max,shape=shape, rate=rate)
	gen.u <- runif(n=1, min=u.min, max=u.max)

    parm$tau <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
	#parm$tau <- 1/sqrt(rgamma(1,shape=shape, rate=rate))
    
    
	#overwrite to avoid zeros and Inf
	if (round(u.min, digits = 5) == 1) # really close to 1
		{parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$min)
		}
	if (round(u.max, digits = 5) == 0) # really close to 0
		{parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$max)
		}
	
	
	parm
	}

########################################

########################################



fn.gen.tau_0  <- function(data, parm)
	{
	###################
	# update tau_0
	###################

    #if(!parm$discrete){
	
	sum.resid.sq <- 0
	count <- 0

	for (g in 1:parm$clust$G)
		{flag.v <- parm$clust$c.v == g
		z.g.v <- parm$clust$s.mt[,g] > 0

		if ((sum(1-z.g.v) > 0) & (sum(flag.v)>0))
			{X.g.mt <- parm$X[!z.g.v,flag.v]
			resid.g.mt <- X.g.mt
			sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
			count <- count + sum(1-z.g.v)*sum(flag.v)
			}
		}

	shape <- 1 + count/2
	rate <- 1 + sum.resid.sq/2

	# minimum possible value of parm$tau_0 = 1.5 * maximum possible value of parm$tau
	# maximum possible value of parm$tau_0 = 3 * sd(as.vector(data$X))
	u.min <- pgamma(1/9 / var(as.vector(data$X),na.rm=TRUE),shape=shape, rate=rate)
	u.max <- pgamma(1/1.5^2/parm$prior$tau.sq$min,shape=shape, rate=rate)
	gen.u <- runif(n=1, min=u.min, max=u.max)

      parm$tau_0 <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))

    #}else{
	
	parm$tau_0 <- 1

	#}
	parm
	}

###############################################################################
# Updating Sd for each covariates 
###############################################################################

fn.disc.tau <- function(parm,data)
	{
	parm$clust$sd.v <- NULL
	
	# data$OX
  	n <- dim(data$OX)[1]
  	p <- dim(data$OX)[2]
	sd  <- matrix(NA,n,p)

	for (xx in 1:p ){
	s.v      <- parm$clust$s.mt[,parm$clust$c.v[xx]]
	theta.v  <- parm$clust$phi.v[s.v]
	
  	# Numerical Stability
  	log.p      <-  -log( 1 + exp(theta.v))  +  theta.v
  	log.var    <-	-2*log(1 + exp(theta.v))  +  theta.v

	sd[,xx]    <- sqrt(1/exp(log.var))
	
	}

	parm$clust$sd.v <- sd
	
parm
}




fn.hyperparameters <- function(data, parm)
	{

	# also updates update tau_int
	parm <- fn.gen.tau(data, parm)
	#parm <- fn.sample.M(parm)
	parm <- fn.DP.hyperparameters(parm,data)

	#parm <- fn.gen.tau_0(data, parm)

	parm

	}


fn.funky <- function(s,t)
	{# on log scale
	lgamma(s+t) - lgamma(s)
	}

fn.d <- function(d, parm)
	{

	# formula in review paper by Lijoi and Prunster
	log.lik <- sum(log(parm$b1 + (1:(parm$clust$G-1))*d)) - fn.funky((parm$b1+1), (parm$p-1)) + sum(fn.funky((1-d), (parm$clust$C.m.vec-1)))

	log.lik
	}




fn.poissonDP.hyperparm <- function(data, parm, w=.01, max.d)
	{


	## update parm$d conditional on parm$b1
	## 1/w must be an integer
	
	
	d.v <- seq(0,max.d,by=w)
	len <- length(d.v)
	d.v <- d.v[-len]
	len <- len-1

	log.lik.v <- sapply(d.v, fn.d, parm)
	# putting 1/2 prior mass on 0 and remaining spread uniformly on positive points in d.v
	log.p.v <- log(.5) + c(0,  rep(-log(len-1),(len-1)))

	
	log.post.v <- log.lik.v + log.p.v
	log.post.v <- log.post.v - max(log.post.v)
	post.v <- exp(log.post.v)
	post.v <- post.v/sum(post.v)

	prop.d <- sample(d.v, size=1, prob=post.v)

	if (prop.d > 0)
		{prop.d <- runif(n=1, min=(prop.d-w), max=(prop.d+w))
		}

	if (prop.d != parm$d)
		{
		# MH ratio for independent proposals and
		# prior same for all d (which is true if 0 wp .5 and \in (0,max.d) wp .5)

		log.ratio <- fn.d(d=prop.d, parm) - fn.d(d=parm$d, parm)
		prob <- min(1, exp(log.ratio))
		flip <- rbinom(n=1, size=1, prob=prob)
		if (flip==1)
			{parm$d <- prop.d
			}
		}

		
	parm

	}


fn.DP.hyperparameters <- function(parm, data){

	post.mu.prec    <- parm$clust$K/parm$clust$tau2^2 + 1/parm$prior$mu2$sd^2
	post.mu.mean	<- 1/post.mu.prec * ( sum(parm$clust$phi.v)/parm$clust$tau2^2 + parm$prior$mu2$mean/parm$prior$mu2$sd^2 )
	parm$clust$mu2  <- rnorm(1,post.mu.mean,sqrt(1/post.mu.prec) )
	#parm$clust$mu2 <- 0
	
	#Updating tau2
	
	post.sigma.shape <- parm$prior$tau2$alpha + parm$clust$K/2
	post.sigma.rate  <- parm$prior$tau2$beta + t(parm$clust$phi.v -parm$clust$mu2)%*%(parm$clust$phi.v -parm$clust$mu2)/2 
	parm$clust$tau2  <- sqrt(1/rgamma(1,post.sigma.shape, scale = 1/post.sigma.rate))
	#parm$clust$tau2  <- 1

parm

}	

fn.log.lik <- function(mu,k, ind, parm){

 ans <- sum( dnorm( parm$X[k,ind], mean = mu, sd = parm$clust$sdm.mt[k,ind] , log =TRUE ) )

 ans
 }

 fn.hyper <- function(k, ind, parm, data){

 post.mean <- sum( data$X[k,ind]/parm$clust$sdm.mt[k,ind]^2 )
 post.prec <- sum( 1/parm$clust$sdm.mt[k,ind]^2 )
 
 out <- list()
 out[[1]] <- post.mean
 out[[2]] <- post.prec
 
 out
 }

########################################

fn.iter <- function(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm,
                    computeMode)
	{

	parm0 <- parm
	parm <- fn.transform(parm,data)
	parm <- fn.offset(parm)    
	
	data$X <- parm$X
		
	parm <- fn.sample.phi(parm,data)

	parm <- fn.equivsd(parm,data)	
	
	parm <- fn.matrix.sd(parm,data)
	
	parm <- fast_PDP_fn.main(parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size, computeMode)

	parm <- fn.element.DP(data, parm, max.row.nbhd.size, row.frac.probes, computeMode)
	
	parm <- fn.poissonDP.hyperparm(data, parm, w=.01, max.d=1)
	
	parm <- fn.hyperparameters(data, parm)
		
	flip <- rbinom(n=1, size=1, prob=.1)
	if (flip==1)
	  {parm <- fn.gen.missing.X(data, parm)
	}

	if (parm$flip.sign)
	  {
	    ############
	    ## Update signs for updated columns
	    ############
	    parm <- PDP_fn.orientation(parm, cc_subset=1:parm$p)
	  }

	if (parm$standardize.X)
	  {parm <- fn.standardize_orient.X(parm)
	  }

	parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)
	if (parm$tBB_flag)
	  {parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt
	  }

	err <- fn.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC: err=",err))
		}

	parm


	}


#####################################################################################################

fn.postBurnin <- function(text, All.Stuff, offset, n.reps, data, parm, dahl.flag, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm)
  {
  
  for (cc in 1:n.reps)
  {
  parm <- fn.iter(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, computeMode)
  All.Stuff$G.v[cc+offset] <- parm$clust$G
  All.Stuff$K.v[cc+offset] <- parm$clust$K
  All.Stuff$tau.v[cc+offset] <- parm$tau
  All.Stuff$tau_0.v[cc+offset] <- parm$tau_0
  All.Stuff$tau_int.v[cc+offset] <- parm$tau_int
  All.Stuff$offset.m[cc+offset,] <- parm$offset
  All.Stuff$mu2.v[cc + offset]   <- parm$clust$mu2
  All.Stuff$tau2.v[cc + offset]  <- parm$clust$tau2

  All.Stuff$tau2.v[cc+offset]   <- parm$clust$tau2
	

  All.Stuff$d.v[cc+offset]   <- parm$d
  All.Stuff$phi.v[cc+offset] <- parm$phi


# summarizing elementwise DP in "fn.groupwise.updates"

  #All.Stuff$row.flip.v[cc+offset]  <- parm$clust$row.flip

  All.Stuff$col_new_clust.v[cc+offset]  <- parm$clust$col.new.flag
  All.Stuff$col_flip.v[cc+offset]  <- parm$clust$col.mh.flip
  All.Stuff$col_exit.v[cc+offset]  <- parm$clust$col.mh.exit

  All.Stuff$nbhd_max[cc+offset] <- round(parm$clust$nbhd_max_dist, digits = 2)

  # least squares calculation
  parm = fn.Dahl(parm, All.Stuff, dahl.flag)
  All.Stuff$mean.taxicab.v[cc+offset] <- mean(true_parm$clust$nbhd.matrix != parm$dahlDist.mt)
  if (!is.finite(All.Stuff$mean.taxicab.v[cc+offset]))
  {stop("Nan")}
  
  if (!dahl.flag)
    {All.Stuff$meanDahlDist.mt = All.Stuff$meanDahlDist.mt + parm$dahlDist.mt
  }
  
  if (dahl.flag)
  {All.Stuff$runningMinDahlDist.v[cc]  <- parm$min_Dahl_dist
  All.Stuff$dahlDist.v[cc]  <- parm$Dahl_dist
  }

  if (cc %% 10 == 0)
    {print(paste(text, "REPS = ",cc,date(),parm$d,All.Stuff$mean.taxicab.v[cc+offset], parm$clust$c.v[6],parm$clust$c.v[47],"***********"))
    }
  
  } # END FOR LOOP
  
  if (!dahl.flag)
  {All.Stuff$meanDahlDist.mt = All.Stuff$meanDahlDist.mt/n.reps
  }
  
  if (dahl.flag)
  {All.Stuff$est_c.v  <- parm$est_c.v
  All.Stuff$min_Dahl_dist = parm$min_Dahl_dist
  }
  
  tmp = list(parm, All.Stuff)
  names(tmp) = c("parm", "All.Stuff")
  
  tmp
  }

##############################################################################

fn.mcmc <- function(text, true, data, n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, dahl.flag=FALSE,
                    standardize.X=FALSE, flip.sign =FALSE, tBB_flag=FALSE, computeMode ="R")
	{

	# initialize
	parm <- fn.init(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm, tBB_flag, standardize.X, flip.sign, computeMode)
	init.parm <- parm
	err <- fn.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC at fn.init: err=",err))
		}

	text="BURNIN..."
	
	for (cc in 1:n.burn)
		{
		parm <- fn.iter(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, computeMode)

		if (cc %% 10 == 0)
			{print(paste(text, "REPS = ",cc,date(),parm$d,"***********"))
			}
		}

	##########################################
	## DEFINE OBJECTS
	##########################################

	All.Stuff <- NULL
	#
	All.Stuff$phi.v <- All.Stuff$d.v <- All.Stuff$tau_0.v <- All.Stuff$tau.v <- All.Stuff$tau_int.v <- All.Stuff$G.v <- All.Stuff$K.v <- array(,2*n.reps)
	All.Stuff$row.flip.v  <- array(,2*n.reps)
	All.Stuff$nbhd_max <- All.Stuff$col_new_clust.v  <- All.Stuff$col_exit.v <- All.Stuff$col_flip.v  <- array(,2*n.reps)

	All.Stuff$mean.taxicab.v  <- array(,2*n.reps)
	All.Stuff$offset.m        <- array(, dim=c(2*n.reps,parm$p ) )
	All.Stuff$mu2.v           <- array(, 2*n.reps)
	All.Stuff$tau2.v           <- array(, 2*n.reps)

	All.Stuff$mean.taxicab.v  <- array(,2*n.reps)
	
	# this is only for n.reps (runs 1 or 2)
	All.Stuff$runningMinDahlDist.v  <- All.Stuff$dahlDist.v  <- array(,n.reps)

	All.Stuff$meanDahlDist.mt <- array(0,c(parm$p,parm$p))
	
	##########################################
	## POST-MCMC RUN 1
	##########################################
	
	text="POST--BURNIN RUN 1..."
	
	tmp = fn.postBurnin(text, All.Stuff, offset=0, n.reps, data, parm, dahl.flag=FALSE, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm)
  	parm = tmp$parm
  	All.Stuff = tmp$All.Stuff

	


	
	##########################################
	## POST-MCMC RUN 2: ALSO COMPUTE LEAST SQUARES ALLOCATION
	##########################################
	
	text="RUN 2: DAHL..."
	
	parm$min_Dahl_dist = parm$p^2
	
	tmp = fn.postBurnin(text, All.Stuff, offset=n.reps, n.reps, data, parm, dahl.flag=TRUE, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm)
	parm = tmp$parm
	All.Stuff = tmp$All.Stuff
	
	#######################
	
	All.Stuff$parm <- parm
	All.Stuff$init.parm <- init.parm

	All.Stuff
	}