
gen.X <- function(n, p, prop.data.type, prop.X.miss, true_parm ) {

	# 6/5/2016 Right now, not concerned about proportion of misses in mixed case
	# prop.X.cont is introduced to allow for a certain proportion of covariates to be continuous
	#X.cont.ind <- rbinom(p, 1 , prop.X.cont)
	X.cont.ind  <- rmultinom(p,1,prop.data.type)
	data.type   <- c(seq(1,5,1))
	X.data.type <- data.type%*%X.cont.ind

	data <- NULL
	data$OX <- data$X <- array(,c(n,p))

	
	for (xx in 1:p)
	{s.v	 <- true_parm$clust$s.mt[,true_parm$clust$c.v[xx]]
	x.v 	 <- theta.v <- array(,n)
	theta.v  <- true_parm$clust$phi.v[s.v]	
		
	if(X.data.type[xx]== 1 ){
		
		#theta.v <- (theta.v - true_parm$clust$mu2)/true_parm$clust$tau2  
		#theta.v <- theta.v +2*(-1)^xx
		x.v     <- pnorm(theta.v) 	
		x.v     <- rbinom(n, 1, x.v)
		#x.v	<- 1/( 1 + exp(-theta.v) ) 	
		#x.v	<- rbinom(n,1,x.v)	 	
	 
	}	
	if(X.data.type[xx]== 7 ){

		theta.v <- theta.v
		x.v    <- rnorm(n, mean = theta.v, sd = true_parm$tau)
		x.v    <- pnorm(x.v)	
		x.v    <- rbinom(n, 1, x.v)	 
	}
	if(X.data.type[xx]== 8 ){
	
		theta.v <- theta.v
		theta_d.v <-true_parm$clust$phi_d.v[true_parm$clust$s.v]
		x.v    <- pnorm(theta_d.v)	
		x.v    <- rbinom(n, 1, x.v)	 
	}		
	
	if(X.data.type[xx]== 2 ){
		
		theta.v <- theta.v + 6.842845
		x.v       <- rnorm(n, mean = theta.v, sd = true_parm$tau)
	}

	if(X.data.type[xx]== 3 ){
		
		 N0      <- 5
		 x.v     <- N0*exp(theta.v)
		 x.v     <- rpois(n , x.v)
		 
	}	
	
	if( X.data.type[xx] == 4){
		

		#theta.v <- theta.v + 1	
		# 5 categories
		true.cutoff <- c(-Inf, -2 , -1, 0, 1, Inf )
		q  <- length(true.cutoff)
		true.cutoff <- matrix( rep(true.cutoff, length(theta.v) ), length(theta.v), q , byrow=TRUE)
		true.prob        <- matrix(NA, length(theta.v), q-1)
		
		for(i in 1:5){
			true.prob[,i] <- pnorm(true.cutoff[,i+1] - theta.v) - pnorm(true.cutoff[,i] - theta.v)
		}
		
		x.v <- theta.v
		for(j in 1: length(theta.v) ){
		x.v[j] <-	which( rmultinom(1, 1 ,true.prob[j,]) == 1)
		}

	
	}
	

		
	#Proportion
	
	if(X.data.type[xx] == 5 ){
		 
		 theta.v <- theta.v
		 phi    <- 24.5
		 mu     <-exp(theta.v)/( 1 + exp(theta.v) )

		 x.v	<- rbeta(n, mu*phi, (1-mu)*phi ) 

		 x.v1   <- x.v
	
		 ind0    <- which(x.v < 0.005 )
		 ind1    <- which(x.v  > 0.995 ) 
		 
		 #for( i in 1:length(ind0) ){
		 #x.v[ind0[i] ] =0.005	
		 #}	
			
		 #for(i in 1:length(ind1) ){
		 #x.v[ind1[i] ] =0.995	
		 #}		

	}	
	
	if(X.data.type[xx] == 6){  
		x.v       <- rnorm(n, mean = theta.v, sd = true_parm$tau)
		x.v       <- exp(x.v)/(1+ exp(x.v) )	
	} 
		
		       
		data$X[,xx] <- x.v
	 	
	}
			
	# Saving data in OX
	data$OX <- data$X
	data$data.type <- X.data.type
    
	###########################################
	# missing X values
	###########################################

	data$num.X.miss <- round(prop.X.miss*n*p)

	if (data$num.X.miss>0)
	  {data$X.missing.x <- sample(1:n,size=data$num.X.miss, replace=TRUE)
	  data$X.missing.y <- sample(1:p,size=data$num.X.miss, replace=TRUE)
	}

	# works even if some (data$X.missing.x,data$X.missing.y) are tied by chance

	for (cc in 1:data$num.X.miss)
	  {data$X[data$X.missing.x[cc],data$X.missing.y[cc]] <- NA
	  }

	###########################################
	# random split of 100 X prop % missing
	###########################################

	n2 <- n
	n1 <- 0
	data$missing.indx <- NULL
	data$non.missing.indx <- 1:n

	data$K.max <- round(n2/2)
	data$G.max <- round(p/2)

	###########################################
	# dummy responses
	###########################################

	data$Y <- rep(0,n2)
	data$delta <- rep(0,n2)
	data$true <- NULL
	data$true$Y <- data$Y
	data$true$delta <- data$delta

	############

	true <- NULL
	true$a.R <- true_parm$clust$M
	true$b0 <- 2.2
	true$b1 <- true_parm$b1

	#########################################
	# generating the R- and C- clusters
	########################################

	true$shift <- 1e-4

	return(list(data = data, true = true))
}
