
fn1.update.element.objects <- function(parm, computeMode)
	{

	parm$clust$A.mt <- parm$clust$s.mt <- array(parm$clust$s.v, c(parm$n2, parm$clust$G))

	for (g in 1:parm$clust$G)
		{
		
		s.g.v <- parm$clust$s.mt[,g]
		s.pos.indx <- s.g.v > 0
		
		if (sum(s.pos.indx) > 0)
			{parm$clust$A.mt[s.pos.indx,g] <- parm$clust$phi.v[s.g.v[s.pos.indx]]
			}
		if ((parm$n2-sum(s.pos.indx)) > 0)
			{parm$clust$A.mt[!s.pos.indx,g] <- 0
			}
		#parm$clust$A.mt[,g] <- parm$clust$phi.v[s.g.v]
		parm$clust$B.mt[,(g+1)] <- parm$clust$A.mt[,g]
		}

	parm$clust$theta.v <- as.vector(parm$clust$A.mt)

	if (computeMode$computeR) {

	  parm$clust$n.vec <- array(,parm$clust$K)
	  for (s in 1:parm$clust$K)
		  {parm$clust$n.vec[s] <- sum(parm$clust$s.v == s)
	  }
	  parm$clust$n0 <- sum(parm$clust$s.v == 0)
	}

	if (computeMode$computeC) {

	  all.n.vec <- .fastTabulateVector(parm$clust$s.v, parm$clust$K, TRUE)

	  if (computeMode$computeR) {
      assertEqual(parm$clust$n0, all.n.vec[1])
	    assertEqual(parm$clust$n.vec, all.n.vec[-1])
	  }

	  parm$clust$n0 <- all.n.vec[1]
	  parm$clust$n.vec <- all.n.vec[-1]
	}

	parm
	}


fn2.update.element.objects <- function(parm, computeMode)
	{

	parm$Y <- parm$X.sd <- array(,c(parm$n2, parm$clust$G))
	parm$Y1 <- parm$X1.sd <- array(NA,c(parm$n2,parm$clust$G))
	parm$Y2 <- parm$X2.sd <- array(NA,c(parm$n2,parm$clust$G))
    parm$Y3 <- parm$X3.sd <- array(NA,c(parm$n2,parm$clust$G))
    parm$Y4 <- parm$X4.sd <- array(NA,c(parm$n2,parm$clust$G))
	parm$Y5 <- parm$X5.sd <- array(NA,c(parm$n2,parm$clust$G))
	parm$Y6 <- parm$X6.sd <- array(NA,c(parm$n2,parm$clust$G))
	
	# group covariate tells which parm$clust$rho.g
	# to use for likelihood calculation
	parm$g <- rep(1:parm$clust$G,each = parm$n2)

	for (g in 1:parm$clust$G)
		{I.g <- (parm$clust$c.v==g)
		 m.g <- parm$clust$C.m.vec[g]

		x.g.v <- x.tmp <- parm$X[,I.g]
		x2.g.v <- x.g.v^2
		if (m.g > 1)
			{x.g.v <- rowMeans(x.tmp)
			x2.g.v <- rowMeans(x.tmp^2)
			}
		 parm$Y[,g] <- x.g.v

		sd.g.v <- rep(0, parm$n2)
		 if (m.g > 1)
			{
			# To make it numerically stable
			err <- 10^-10
			sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1) + err )
			}
		parm$X.sd[,g] <- sd.g.v
		
		####################################################
		### Split nbhd
		####################################################
		
		### Discrete
		### Data Type
		
		#Probit
		I1.g <-(parm$clust$c.v == g & parm$data.type == 1)
		
		if(sum(I1.g) > 0 ){
		
		x.g.v <- x.tmp <- parm$X
		x.g.v <- x.tmp <- parm$X[,I1.g]
		
		x2.g.v <- x.g.v^2
		m.g <- sum(I1.g)
		if (m.g > 1 )
		{
			x.g.v <- rowMeans(x.tmp)
			x2.g.v <- rowMeans(x.tmp^2)
		}
		parm$Y1[,g] <- x.g.v
		
		sd.g.v <- rep(0, parm$n2)
		 if (m.g > 1)
		 	{# To make it numerically stable
			 err <- 10^-10
			sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1) + err)
			}
		parm$X1.sd[,g] <- sd.g.v
		}
		
		### Continuous
		I2.g <-(parm$clust$c.v == g & parm$data.type == 2)
		if( sum(I2.g) > 0 ) {
		x.g.v <- x.tmp <- parm$X[,I2.g]
		x2.g.v <- x.g.v^2
		m.g <- sum(I2.g)
		if (m.g > 1 )
		{
			x.g.v <- rowMeans(x.tmp)
			x2.g.v <- rowMeans(x.tmp^2)
		}
		parm$Y2[,g] <- x.g.v
		
		sd.g.v <- rep(0, parm$n2)
		 if (m.g > 1)
			{
			# To make it numerically stable
			err <- 10^-10
			sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1) + err)
			}
		parm$X2.sd[,g] <- sd.g.v
		}

		### Poisson
		I3.g <-(parm$clust$c.v == g & parm$data.type == 3)
		if( sum(I3.g) > 0 ) {
		x.g.v <- x.tmp <- parm$X[,I3.g]
		x2.g.v <- x.g.v^2
		m.g <- sum(I3.g)
		if (m.g > 1 )
		{
			x.g.v <- rowMeans(x.tmp)
			x2.g.v <- rowMeans(x.tmp^2)
		}
		parm$Y3[,g] <- x.g.v
		
		sd.g.v <- rep(0, parm$n2)
		 if (m.g > 1)
			{
			# To make it numerically stable
			err <- 10^-10
			sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1) + err)
			}
		parm$X3.sd[,g] <- sd.g.v
		}		
		
		### Ordinal
		I4.g <-(parm$clust$c.v == g & parm$data.type == 4)
		if( sum(I4.g) > 0 ) {
		x.g.v <- x.tmp <- parm$X[,I4.g]
		x2.g.v <- x.g.v^2
		m.g <- sum(I4.g)
		if (m.g > 1 )
		{
			x.g.v <- rowMeans(x.tmp)
			x2.g.v <- rowMeans(x.tmp^2)
		}
		parm$Y4[,g] <- x.g.v
		
		sd.g.v <- rep(0, parm$n2)
		 if (m.g > 1)
			{
			# To make it numerically stable
			err <- 10^-10
			sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1) + err)
			}
		parm$X4.sd[,g] <- sd.g.v
		}
		
		### Proportion
		I5.g <-(parm$clust$c.v == g & parm$data.type == 5)
		if( sum(I5.g) > 0 ) {
		x.g.v <- x.tmp <- parm$X[,I5.g]
		x2.g.v <- x.g.v^2
		m.g <- sum(I5.g)
		if (m.g > 1 )
		{
			x.g.v <- rowMeans(x.tmp)
			x2.g.v <- rowMeans(x.tmp^2)
		}
		parm$Y5[,g] <- x.g.v
		
		sd.g.v <- rep(0, parm$n2)
		 if (m.g > 1)
			{
			# To make it numerically stable
			err <- 10^-10
			sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1) + err)
			}
		parm$X5.sd[,g] <- sd.g.v
		}
		
		### Continuous Proportion
		I6.g <-(parm$clust$c.v == g & parm$data.type == 6)
		if( sum(I6.g) > 0 ) {
		x.g.v <- x.tmp <- parm$X[,I6.g]
		x2.g.v <- x.g.v^2
		m.g <- sum(I6.g)
		if (m.g > 1 )
		{
			x.g.v <- rowMeans(x.tmp)
			x2.g.v <- rowMeans(x.tmp^2)
		}
		parm$Y6[,g] <- x.g.v
		
		sd.g.v <- rep(0, parm$n2)
		 if (m.g > 1)
			{
			# To make it numerically stable
			err <- 10^-10
			sd.g.v <- sqrt((x2.g.v - x.g.v^2)*m.g/(m.g-1) + err)
			}
		parm$X6.sd[,g] <- sd.g.v
		}
		
		}
		

	parm$N <- parm$n2 * parm$clust$G

	parm$Y <- as.vector(parm$Y)

	parm$X.sd <- as.vector(parm$X.sd)

	parm$Y1 <- as.vector(parm$Y1)
	
	parm$X1.sd <- as.vector(parm$X1.sd)
	
	parm$Y2  <- as.vector(parm$Y2)
	
	parm$X2.sd <- as.vector(parm$X2.sd)
	
	
	parm$Y3  <- as.vector(parm$Y3)
	
	parm$X3.sd <- as.vector(parm$X3.sd)
	
	
	parm$Y4  <- as.vector(parm$Y4)
	
	parm$X4.sd <- as.vector(parm$X4.sd)
	
	
	parm$Y5  <- as.vector(parm$Y5)
	
	parm$X5.sd <- as.vector(parm$X5.sd)

	parm$Y6  <- as.vector(parm$Y6)
	
	parm$X6.sd <- as.vector(parm$X6.sd)
	#####################

	parm <- fn1.update.element.objects(parm, computeMode)

	parm
	}



fn.element.DP <- function(data, parm, max.row.nbhd.size, row.frac.probes,
                          computeMode)
{

  if (parm$standardize.X)
  {parm <- fn.standardize_orient.X(parm)
  }

	# essentially, a Bush-Mac move: given groups, the parm$N=n2XG number of
	# invidividual elements (summaries of microarray elements) belonging to group g>0
	# are updated for s (phi) and z

	#print(7)
	parm <- fn2.update.element.objects(parm, computeMode)
	
	#Update sd
	parm <- fn.equivsd(parm,data)
	parm <- fn.matrix.sd(parm,data)
	
	#print(8)
	#print(mean(data$OX))
	parm <- element_fn.fast.DP(parm, data, max.row.nbhd.size, row.frac.probes, computeMode)

	#print(9)
	#############################
	## Important: do not remove call to fn1.update.element.objects
	## updates A.mt, theta.v, B.mt, tBB.mt, s.v, s.mt, n.vec, n0
	#############################

	parm <- fn1.update.element.objects(parm, computeMode)

  	parm
}



