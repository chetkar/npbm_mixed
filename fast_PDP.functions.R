
PDP_fn.compact.consistency.check <- function(parm)
	{err <- 0

	if (max(sort(unique(parm$clust$c.v))) != parm$clust$G)
		{err <- 1
		}

	if (sum(parm$clust$C.m.vec==0) > 0)
		{err <- 1.5
		}

	err <- PDP_fn.consistency.check(parm)

	err
	}

PDP_fn.check.nbhd <- function(parm)
{
  max_dist.v <- array(,length(parm$clust$col.nbhd.k))

  for (zz in 1:length(parm$clust$col.nbhd.k))
  {
    v1 <- parm$clust$col_post.prob.mt[,parm$clust$col.nbhd.k[zz]]
    m1 <- matrix(parm$clust$col_post.prob.mt[,parm$clust$col.nbhd[[zz]]],ncol=length(parm$clust$col.nbhd[[zz]]))
    max_dist.v[zz] <- max(2*(1-colSums(sqrt(v1*m1))))
  }

  parm$clust$nbhd_max_dist <- max(max_dist.v)

  parm

}

PDP_fn.consistency.check <- function(parm)	{

err <- 0


	#if (sum(parm$clust$C.m.vec) + parm$clust$C.m0 != parm$p)
	#	{err <- 4
	#	}

	if (length(parm$clust$n.vec) != parm$clust$K)
		{err <- 9
		}

	if (length(parm$clust$phi.v) != parm$clust$K)
		{err <- 10
		}
	
	err
}


PDP_fn.swap.clusters <- function(parm, g1, g2, computeMode)
	{

	####################################################
	# swap the group labels g1 and g2
	####################################################

	ind1 <- parm$clust$c.v == g1
	ind2 <- parm$clust$c.v == g2
	parm$clust$c.v[ind1] <- g2
	parm$clust$c.v[ind2] <- g1
	
	if (computeMode$computeR) {
	  buffer <- parm$clust$s.mt[,g1]
	  parm$clust$s.mt[,g1] <- parm$clust$s.mt[,g2] # HOT 3
	  parm$clust$s.mt[,g2] <- buffer
	} else {
	  .swapIntegerMatrix(parm$clust$s.mt, g1, g2, FALSE)
	}

	buffer <- parm$clust$beta.v[g1]
	parm$clust$beta.v[g1] <- parm$clust$beta.v[g2]
	parm$clust$beta.v[g2] <- buffer

	buffer <- parm$clust$gamma.v[g1]
	parm$clust$gamma.v[g1] <- parm$clust$gamma.v[g2]
	parm$clust$gamma.v[g2] <- buffer

	buffer <- parm$clust$C.m.vec[g1]
      parm$clust$C.m.vec[g1] <- parm$clust$C.m.vec[g2]
      parm$clust$C.m.vec[g2] <- buffer

	buffer <- parm$clust$small.indx[g1]
      parm$clust$small.indx[g1] <- parm$clust$small.indx[g2]
      parm$clust$small.indx[g2] <- buffer

	buffer <- parm$clust$order.v[g1]
      parm$clust$order.v[g1] <- parm$clust$order.v[g2]
      parm$clust$order.v[g2] <- buffer

	#####################

	buffer <- parm$clust$A.mt[,g1]
	parm$clust$A.mt[,g1] <- parm$clust$A.mt[,g2]
	parm$clust$A.mt[,g2] <- buffer

	#buffer <- parm$clust$B.mt[,(g1+1)]
	#parm$clust$B.mt[,(g1+1)] <- parm$clust$B.mt[,(g2+1)]
	#parm$clust$B.mt[,(g2+1)] <- buffer

	if (parm$tBB_flag)
	{
	if (computeMode$computeR) {

	  # first swap columns
	  buffer <- parm$clust$tBB.mt[,(g1+1)]
	  parm$clust$tBB.mt[,(g1+1)] <- parm$clust$tBB.mt[,(g2+1)] # HOT 3
	  parm$clust$tBB.mt[,(g2+1)] <- buffer
	  # then swap rows
	  buffer <- parm$clust$tBB.mt[(g1+1),]
	  parm$clust$tBB.mt[(g1+1),] <- parm$clust$tBB.mt[(g2+1),]
	  parm$clust$tBB.mt[(g2+1),] <- buffer

	} else {
	  .swap(parm$clust$tBB.mt, g1 + 1, g2 + 2, TRUE)
	}
	} # end if (parm$tBB_flag)

	parm
	}

PDP_fn.clip.clusters <- function(parm, keep)
	{

	parm$clust$s.mt <- parm$clust$s.mt[,keep]

	parm$clust$beta.v <- parm$clust$beta.v[keep]
	parm$clust$gamma.v <- parm$clust$gamma.v[keep]
	parm$clust$C.m.vec <- parm$clust$C.m.vec[keep]
      parm$clust$small.indx <- parm$clust$small.indx[keep]
      parm$clust$order.v <- parm$clust$order.v[keep]

     	parm$clust$A.mt <- parm$clust$A.mt[,keep]

	indx2 <- c(1,(keep+1))
	parm$clust$B.mt <- parm$clust$B.mt[,indx2]

	if(parm$tBB_flag)
	{parm$clust$tBB.mt <- parm$clust$tBB.mt[indx2,indx2]
	}

	parm
	}



###########################################################


PDP_fn.log.lik <- function(gg, col, parm, data,colSums=TRUE )
	{

	x.mt <- parm$X[,col]
	
	# This might be better but it is computationally costly.
	
	a2.v <- parm$clust$A.mt[,gg]
	z.g.v <- parm$clust$s.mt[,gg] > 0
	log.lik.v <- rep(0,length(col)) ## HOT

	a2.1.v <- parm$clust$A.mt[,gg]
	small.X.1 <- matrix(x.mt, ncol=length(col)) # HOT
	
	#parm <- fn.equivsd(parm,data)
	sd    <- parm$clust$sdm.mt[,col]
		
	tmp <- colSums(-.5*(small.X.1 - a2.1.v)^2/sd^2)
	log.lik.v   <-  tmp 
	log.lik.v <- log.lik.v - parm$n2*(.5*log(2*pi) + mean(log(sd)) )
	
	log.lik.v
	}

PDP_fn.nbhd <- function(relative_I, parm, max.col.nbhd.size)
{
  if (length(relative_I)>1)
  {relative_k <- sample(relative_I, size=1)
  }
  if (length(relative_I)==1)
  {relative_k <- relative_I
  }

  post.prob.mt <- parm$clust$col_subset_post.prob.mt

  tmp1.mt <- matrix(post.prob.mt[,relative_I], ncol = length(relative_I)) ## HOT
  tmp2.v <- post.prob.mt[,relative_k]
  tmp3.mt <- sqrt(tmp1.mt * tmp2.v) ## HOT
  H.v <-  2 * (1 - colSums(tmp3.mt))

  cutoff <- parm$col.delta
  flag.v <- which(H.v <= cutoff)
  relative_I.k <- relative_I[flag.v]

  if (length(relative_I.k) > max.col.nbhd.size) {
    relative_I.k <- relative_I[rank(H.v, ties="random") <= max.col.nbhd.size]
  }

  relative_I.k <- sort(relative_I.k)

  relative_I <- sort(setdiff(relative_I, relative_I.k))
  relative_I <- sort(relative_I)

  list(relative_k, relative_I.k, relative_I)

}


PDP_fn.post.prob.and.delta <- function(parm,data, max.col.nbhd.size, col.frac.probes, computeMode) {

  col.subset <- 1:parm$p

  if (computeMode$computeR) {

    ################################################
    ### Compute pmf of cluster variables w_1,...,w_p
    ###############################################

    #prior.prob.v <- c(parm$clust$C.m0, parm$clust$C.m.vec)
    prior.prob.v <- c(parm$clust$C.m.vec)
	small <- 1e-3 # compared to 1
    prior.prob.v[prior.prob.v < small] <- small

    subset_log.ss.mt <- array(, c((parm$clust$G ), length(col.subset)))

	 for (gg in 1:parm$clust$G)
    {subset_log.ss.mt[(gg),] <- PDP_fn.log.lik(gg, col = col.subset, parm, data )
    }

    subset_log.ss.mt <- subset_log.ss.mt + log(prior.prob.v)

    maxx.v <- apply(subset_log.ss.mt, 2, max)

    subset_log.ss.mt <- t(t(subset_log.ss.mt) - maxx.v)
    subset_ss.mt <- exp(subset_log.ss.mt)

    col.sums.v <- colSums(subset_ss.mt)
    subset_ss.mt <- t(t(subset_ss.mt)/col.sums.v)

    # replace zeros by "small"
    small2 <- 1e-5
    subset_ss.mt[subset_ss.mt < small2] <- small2

    # again normalize
    col.sums.v <- colSums(subset_ss.mt)
    subset_ss.mt <- t(t(subset_ss.mt)/col.sums.v)

    parm$clust$col_post.prob.mt <- array(,c((parm$clust$G), parm$p))
    parm$clust$col_post.prob.mt[,col.subset] <- subset_ss.mt

    dimnames(parm$clust$col_post.prob.mt) <- list(1:parm$clust$G, 1:parm$p)

    parm$clust$col_subset_post.prob.mt <- subset_ss.mt
    dimnames(parm$clust$col_subset_post.prob.mt) <- list(1:parm$clust$G, 1:length(col.subset))

    #########################################
    ### now compute the delta-neighborhoods
    #########################################

    if (computeMode$computeC) { # debugging
      savedSeed <- .GlobalEnv$.Random.seed # For debugging purposed only
    }

    parm$clust$col.nbhd <- NULL
    parm$clust$col.nbhd.k <- NULL
    relative_I <- 1:length(col.subset)

    while (length(relative_I) >= 1) {
      tmp <- PDP_fn.nbhd(relative_I, parm, max.col.nbhd.size) # HOT inside function
      relative_k <- tmp[[1]]
      relative_I.k <- tmp[[2]]
      relative_I <- tmp[[3]]
      #
      parm$clust$col.nbhd <- c(parm$clust$col.nbhd, list(col.subset[relative_I.k]))
      parm$clust$col.nbhd.k <- c(parm$clust$col.nbhd.k, col.subset[relative_k])
    }

    parm <- PDP_fn.check.nbhd(parm)

  }

  if (computeMode$computeC) {

    if (computeMode$computeR) { # debugging
      .GlobalEnv$.Random.seed <- savedSeed # Roll back PRNG
    }

    test <- .computeColumnPmfAndNeighborhoods(computeMode$device$engine,
                                        parm$clust$C.m0, parm$clust$C.m.vec, 1e-3, 1e-5,
                                        parm$clust$G, parm$n2,
                                        parm$Y, parm$X,
                                        parm$clust$A.mt, parm$clust$s.mt,
                                        col.subset,
                                        parm$clust$C.m.vec, parm$p,
                                        parm$clust$phi.v, parm$tau, parm$tau_0, parm$tau_int,
                                        max.col.nbhd.size, parm$col.delta, TRUE)

    if (computeMode$computeR) { # debugging
      assertEqual(test$index, parm$clust$col.nbhd.k)
      assertEqual(test$neighbor, unlist(parm$clust$col.nbhd))
      assertEqual(test$neighborhoodMax, parm$clust$nbhd_max_dist, computeMode$tolerance)
    }

    # Convert from simple flat format to list of int vectors
    end <- test$offset[-1] - 1
    begin <- test$offset
    length(begin) <- length(begin) - 1

    parm$clust$col.nbhd <- lapply(1:length(begin),
                                  FUN = function(x) {
                                    test$neighbor[begin[x]:end[x]]
                                  })
    parm$clust$col.nbhd.k <- test$index
    parm$clust$nbhd_max_dist <- test$neighborhoodMax

  } # computeMode

  ## END

  parm
}


###########################################################

PDP_fn.gibbs <- function(k, parm, data, computeMode)
{	k <- parm$k

	err <- PDP_fn.consistency.check(parm)
	if (err > 0)
		{stop(paste("GIBBS - 0: failed consistency check: err=",err))
		}

	old.c.k <- parm$clust$c.v[k]

	###############

	if (old.c.k > 0)
		{parm$clust$C.m.vec[old.c.k] <- parm$clust$C.m.vec[old.c.k] - 1
		}
	
	x.mt <- matrix(parm$X[,k], ncol=1)


	if (computeMode$computeR) {

	# intercept cluster or any existing cluster
	L.v <- sapply(1:parm$clust$G, PDP_fn.log.lik, k, parm, data)

    }
	

	if (computeMode$computeC) {

    test <- .computePdpLogLikelihood(computeMode$device$engine, k, parm$X,
                                     parm$clust$A.mt, parm$clust$s.mt,
                                     parm$clust$G, parm$n2,
                                     parm$tau, parm$tau_0, parm$tau_int, FALSE)

    if (computeMode$computeR) { # debugging
      assertEqual(test$logLikelihood, L.v, computeMode$tolerance)
    }

    L.v <- test$logLikehood
  } # computeMode
  # NB: returned logLikelihood differ from those computed above by approx 1e-15.  I believe this is due to non-transitivity of FLOPs

	#######################################################
	### emptied clusters are gone forever under Gibbs sampling
	#######################################################

	emptied.indx <- which(parm$clust$C.m.vec==0)
	new.G <- parm$clust$G - length(emptied.indx)

  if (length(emptied.indx) >0)
  	{
    	if (computeMode$computeR) {
    	  new.s.mt <- parm$clust$s.mt[,-emptied.indx] # HOT 3
    	  new.n.vec <- array(,parm$clust$K)
    	  for (pp in 1:parm$clust$K) {
    	    new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
    	  }
     	}

      if (computeMode$computeC) {

    	  tab <- .fastTabulateExcludeEmptiedIndices(parm$clust$s.mt, emptied.indx, parm$clust$K, TRUE)

    	  if (computeMode$computeR) { # debugging
    	    assertEqual(tab[1], new.n0)
          assertEqual(tab[-1], new.n.vec)
    	  }

    	  new.n0 <- tab[1]
    	  new.n.vec <- tab[-1]

    	} # computeMode

    	emptied.s.indx <- which(new.n.vec == 0)
    	new.K <- parm$clust$K - length(emptied.s.indx)

   	}

  if (length(emptied.indx) ==0)
  	{
    	new.s.mt <- parm$clust$s.mt
    	new.n.vec <- parm$clust$n.vec
    	emptied.s.indx <- which(new.n.vec==0)
    	new.K <- parm$clust$K
    #	new.n0 <- parm$clust$n0
    }

  ## generate auxilliary P vector

  tmp.M <- rep(parm$clust$M/new.K,parm$clust$K)
  tmp.alpha <- tmp.M+new.n.vec
  tmp.alpha[emptied.s.indx] <- 0
  P.aux <- rgamma(parm$clust$K,c(tmp.alpha),1)
  P.aux <- P.aux/sum(P.aux)

  ## marginal likelihood of new cluster

  if (computeMode$computeR) {

    marg.log.lik.v <- array(,length(x.mt))
    sd  <- parm$clust$sdm.mt[,k]
	
	for (tt in 1:length(x.mt))
    {
  
		tmp.log.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=sd[tt],log=TRUE) # HOT 3
		tmp.lik.v   <- exp(tmp.log.lik.v)
		tmp.marg.v  <- tmp.lik.v*P.aux
		
		marg.log.lik.v[tt] <- log(sum(tmp.marg.v))
    }
      marg.log.lik <- sum(marg.log.lik.v)

  }

  if (computeMode$computeC) {

    test <- .computeMarginalLikelihood(computeMode$device$engine,
                                       x.mt,
                                       parm$clust$phi.v,
                                       P.aux,
                                       parm$tau, parm$tau_0,
                                       FALSE, # no sampling
                                       computeMode$exactBitStream)

    if (computeMode$computeR) { # debugging
      assertEqual(test$logMarginalLikelihood, marg.log.lik, computeMode$tolerance)
    }

    marg.log.lik <- test$logMarginalLikelihood

  }

  L.v <- c(L.v, marg.log.lik)

	log.prior.v <- array(NA, (1+parm$clust$G))

	if (length(emptied.indx) >0)
		{
		#log.prior.v[-(emptied.indx+1)] <- log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec[-emptied.indx]-parm$d), (parm$b1+new.G*parm$d)))
		log.prior.v[-(emptied.indx)] <- log(c((parm$clust$C.m.vec[-emptied.indx]-parm$d), (parm$b1+new.G*parm$d)))
		log.prior.v[emptied.indx] <- -Inf
		}

	if (length(emptied.indx) ==0)
		{#log.prior.v <- log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec-parm$d), (parm$b1+new.G*parm$d)))
		 log.prior.v <- log(c((parm$clust$C.m.vec-parm$d), (parm$b1+new.G*parm$d)))
		}

	tmp2 <- log.prior.v + L.v
	maxx <- max(tmp2)
	tmp2 <- tmp2 - maxx

	tmp2 <- exp(tmp2)

	parm$clust$post.k <- tmp2
	parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)

	########################################################################
	# store current state
	old.parm <- parm

	########################################################################

	new.c.k <- sample(1:(parm$clust$G+1), size=1, replace=TRUE, prob=parm$clust$post.k)
	parm$clust$c.v[k] <- new.c.k
	new.flag <- new.c.k == (parm$clust$G+1)

	#######################


	if (!new.flag)
		{parm$clust$C.m.vec[new.c.k] <- parm$clust$C.m.vec[new.c.k] + 1
		}

	if (new.flag)
  {
    ###generate the latent vector first, condition on the single kth column
    cand.s.v.k <- array(,length(x.mt))
	sd <- parm$clust$sdm.mt[,k]

    for (tt in 1:length(x.mt))
	  {
	   
			tmp.log.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=sd[tt],log=TRUE) # HOT 3
			tmp.lik.v   <- exp(tmp.log.lik.v)
			tmp.prob.v <- tmp.lik.v*P.aux
			prob.gen.v  <- tmp.prob.v/sum(tmp.prob.v)
			cand.s.v.k[tt]<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=prob.gen.v) # HOT
    } # end for loop

    parm$cand$s.v.k <- cand.s.v.k

    parm$cand$n.vec.k <- array(,parm$clust$K)
    for (gg in 1:parm$clust$K)
    	{parm$cand$n.vec.k[gg] <- sum(cand.s.v.k==gg)
    	}
   
   ##################
   parm$clust$G <- parm$clust$G + 1
   parm$clust$C.m.vec <- c(parm$clust$C.m.vec, 1)
   parm$clust$beta.v <- c(parm$clust$beta.v, 0)
   parm$clust$gamma.v <- c(parm$clust$gamma.v,0)
   parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)
   parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$n2)
   parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k
   parm$N  <- sum(parm$clust$n.vec)
   tmp.a.v <- array(,parm$n2)
   s.G.v <- parm$cand$s.v.k
   tmp.a.v         <- parm$clust$phi.v[s.G.v]
   parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT
   parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v) # HOT
   parm            <- fn.margsd(parm,k)
   
   
		#Reason:
		#Update sd
		#parm <- fn.equivsd(parm,data)
		#parm <- fn.matrix.sd(parm,data)		
		
		if (parm$tBB_flag)
		{
		if (computeMode$computeR) {
		  parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt # HOT
		} else {
		  parm$clust$tBB.mt <- .fastXtX(parm$clust$B.mt)
		}
		}

  } # end  if (new.flag)


	list(parm, new.flag)
}

###########################################################

PDP_fn.fast_col <- function(cc, parm, data, computeMode)
{
  k <- parm$k <- parm$clust$col.nbhd.k[[cc]]
  I.k <- parm$clust$col.nbhd[[cc]]

   err <- PDP_fn.consistency.check(parm)
   if (err > 0)
   {stop(paste("FAST - 0: failed consistency check: err=",err))
   }

  # store so that we can revert to this state if MH propsal is rejected
  init.cc.parm <- parm

  old.c.k <- parm$clust$c.v[I.k]


  if (computeMode$computeR) {

    parm$clust$C.m.vec.k <- array(,parm$clust$G)
    for (gg in 1:parm$clust$G) {
      parm$clust$C.m.vec.k[gg] <- sum(old.c.k == gg)
    }

  }

  if (computeMode$computeC) {

    all.n.vec <- .fastTabulateVector(old.c.k, parm$clust$G, TRUE)

    if (computeMode$computeR) { # debugging
      assertEqual(all.n.vec[1], parm$clust$C.m0.k)
      assertEqual(all.n.vec[-1], parm$clust$C.m.vec.k)
    }

    parm$clust$C.m0.k <- all.n.vec[1] # test1
    parm$clust$C.m.vec.k <- all.n.vec[-1] # test2
  }

  parm$clust$C.m.vec.k.comp <- parm$clust$C.m.vec - parm$clust$C.m.vec.k

  x.mt <- matrix(parm$X[,k], ncol=1)

  if (computeMode$computeR) {

	L.v  <- sapply(1:parm$clust$G, PDP_fn.log.lik, k, parm, data) # HOT
  }

  if (computeMode$computeC) {
    test <- .computePdpLogLikelihood(computeMode$device$engine, k, parm$X,
                                     parm$clust$A.mt, parm$clust$s.mt,
                                     parm$clust$G, parm$n2,
                                     parm$tau, parm$tau_0, parm$tau_int, FALSE)

    if (computeMode$computeR) { # debugging
      assertEqual(test$logLikelihood, L.v, computeMode$tolerance)
    }

    L.v <- test$logLikehood
  } # computeMode
  # NB: returned logLikelihood differ from those computed above by approx 1e-15.  I believe this is due to non-transitivity of FLOPs

   emptied.indx <- which(parm$clust$C.m.vec.k.comp==0)
   new.G <- parm$clust$G - length(emptied.indx)

  if (length(emptied.indx) >0)
  {
    if (computeMode$computeR) {
      new.s.mt <- parm$clust$s.mt[,-emptied.indx] # HOT 3
      new.n.vec <- array(,parm$clust$K)
      for (pp in 1:parm$clust$K) {
        new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
      }
     # new.n0 <- sum(new.s.mt == 0) # HOT 3
    }

    if (computeMode$computeC) { # computeMode
      tab <- .fastTabulateExcludeEmptiedIndices(parm$clust$s.mt, emptied.indx, parm$clust$K, TRUE)

      if (computeMode$computeR) { # debugging
        assertEqual(tab[1], new.n0)
        assertEqual(tab[-1], new.n.vec)
      }

      new.n0 <- tab[1]
      new.n.vec <- tab[-1]

    } # computeMode

    emptied.s.indx <- which(new.n.vec == 0)
    new.K <- parm$clust$K - length(emptied.s.indx)
	
  }

  if (length(emptied.indx) ==0)
  {
    new.s.mt <- parm$clust$s.mt
    new.n.vec <- parm$clust$n.vec
    emptied.s.indx <- which(new.n.vec==0)
    new.K <- parm$clust$K
   # new.n0 <- parm$clust$n0
  }

  ## generate auxilliary P vector

  tmp.M <- rep(parm$clust$M/new.K,parm$clust$K)
  tmp.alpha <- tmp.M+new.n.vec
  tmp.alpha[emptied.s.indx] <- 0
  P.aux <- rgamma(parm$clust$K,c(tmp.alpha),1)
  P.aux <- P.aux/sum(P.aux)

  # START

  if (computeMode$computeR) {

    if (computeMode$computeC) { # debugging
      savedSeed <- .GlobalEnv$.Random.seed # For debugging purposed only
    }

  ## marginal likelihood of new cluster
  marg.log.lik.v <- cand.s.v.k <- array(,length(x.mt))

   sd  <- parm$clust$sdm.mt[,k]
   for (tt in 1:length(x.mt))
	  {
	   
		tmp.log.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=sd[tt],log=TRUE) # HOT 3
		tmp.lik.v   <- exp(tmp.log.lik.v)
		
		tmp.marg.v <- tmp.lik.v*P.aux
		#Numerical Stability
	    tmp.marg.v[tmp.marg.v < 10^-10] <- 10^-10
		marg.log.lik.v[tt] <- log(sum(tmp.marg.v))
        cand.s.v.k[tt]<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v)
  }

  parm$cand$s.v.k <- cand.s.v.k
  marg.log.lik <- sum(marg.log.lik.v)

  # SWTICH
  }

  if (computeMode$computeC) {

    if (computeMode$computeR) { # debugging
      .GlobalEnv$.Random.seed <- savedSeed # Roll back PRNG
    }

  test <- .computeMarginalLikelihood(computeMode$device$engine,
                                     x.mt,
                                     parm$clust$phi.v,
                                     P.aux,
                                     parm$tau, parm$tau_0,
                                     TRUE, # do sampling
                                     computeMode$exactBitStream)

    if (computeMode$computeR) { # debugging
      assertEqual(test$logMarginalLikelihood, marg.log.lik, computeMode$toleranace)
      assertEqual(test$sVk, cand.s.v.k)
    }

  parm$cand$s.v.k <- test$sVk
  marg.log.lik <- test$logMarginalLikelihood

  }

  # END

  L.v <- c(L.v, marg.log.lik)

  ##
   log.prior.v <- array(NA, (1+parm$clust$G))

   spread.mass <- (parm$b1+new.G*parm$d)/(1+length(emptied.indx))
   
   if (length(emptied.indx) >0)
   {
    #log.prior.v[-(emptied.indx+1)] <- log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec[-emptied.indx]-parm$d), spread.mass))
    #log.prior.v[emptied.indx+1] <- log(spread.mass)
    log.prior.v[-(emptied.indx)] <- log(c( (parm$clust$C.m.vec[-emptied.indx]-parm$d), spread.mass))
    log.prior.v[emptied.indx] <- log(spread.mass)
   }

   if (length(emptied.indx) ==0)
   {#log.prior.v <- log(c((parm$b0+parm$clust$C.m0.k.comp), (parm$clust$C.m.vec.k.comp-parm$d), spread.mass))
	log.prior.v <- log(c((parm$clust$C.m.vec.k.comp-parm$d), spread.mass))
   }

   tmp2 <- log.prior.v + L.v
   maxx <- max(tmp2)
   tmp2 <- tmp2 - maxx

   tmp2 <- exp(tmp2)

   parm$clust$post.k <- tmp2
   parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)

   ########################################################################
   # store current state
   old.parm <- parm

   ########################################################################
   #Numerical Stability
   #parm$clust$post.k[parm$clust$post.k < 10^-10] <- 10^-10

   new.c.k <- sample(1:(parm$clust$G+1), size=length(I.k), replace=TRUE, prob=parm$clust$post.k)

  exit <- (sum(new.c.k != old.c.k)==0)
  flip <- TRUE
  new.flag <- FALSE

  if (!exit) # CONTINUE W/O EXITING FUNCTION
  {

    parm$clust$c.v[I.k] <- new.c.k
	new.count <- sum(new.c.k == (parm$clust$G+1))
	
	#######################
	new.prop <- 0

    if (computeMode$computeR) {

      for (gg in 1:parm$clust$G)
      {count.gg <- sum(new.c.k==gg)
       parm$clust$C.m.vec[gg] <- parm$clust$C.m.vec.k.comp[gg] + count.gg

       if (count.gg > 0)
          {new.prop <- new.prop + log(parm$clust$post.k[gg])*count.gg
          }
      }

      }

      if (computeMode$computeC) {

        all.count <- .fastTabulateVector(new.c.k, parm$clust$G + 1, TRUE)
        test1 <- parm$clust$C.m0.k.comp + all.count[1] # test1
        test2 <- parm$clust$C.m.vec.k.comp + all.count[2:(parm$clust$G + 1)] # test2
        test3 <- .fastSumSafeLog(parm$clust$post.k, all.count, parm$clust$G + 1) # test3

        if (computeMode$computeR) { # debugging
          assertEqual(test1, parm$clust$C.m0)
          assertEqual(test2, parm$clust$C.m.vec)
          assertEqual(test3, new.prop)
        }

        parm$clust$C.m0 <- test1
        parm$clust$C.m.vec <- test2
        new.prop <- test3


      }

    if (new.count > 0)
      {
      parm$clust$G <- parm$clust$G + 1
	  parm$clust$C.m.vec <- c(parm$clust$C.m.vec, new.count)
	  new.prop <- new.prop + log(parm$clust$post.k[parm$clust$G])*new.count
	  parm$cand$n.vec.k <- array(,parm$clust$K)
        for (ss in 1:parm$clust$K)
        {parm$cand$n.vec.k[ss] <- sum(parm$cand$s.v.k==ss)
        }
      ##################


        parm$clust$beta.v <- c(parm$clust$beta.v, 0)
        parm$clust$gamma.v <- c(parm$clust$gamma.v,0)

        parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)

        parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$n2)

        parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k
 		parm$N  <- sum(parm$clust$n.vec)

        tmp.a.v <- array(,parm$n2)
        s.G.v <- parm$cand$s.v.k
        #indxx <- s.G.v==0
        #tmp.a.v[indxx] <- 0
        #tmp.a.v[!indxx] <- parm$clust$phi.v[s.G.v[!indxx]] 
		tmp.a.v <- parm$clust$phi.v[s.G.v]
        #
        parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT
        parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v)
		parm            <- fn.margsd(parm,I.k)

		#Update sd
		#parm <- fn.equivsd(parm,data)
		#parm <- fn.matrix.sd(parm,data)
		
        
		if (parm$tBB_flag)
        {
        if (computeMode$computeR) {
          parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt # HOT
        } else {
          parm$clust$tBB.mt <- .fastXtX(parm$clust$B.mt)
        }
        }

      } # end  if (new.count > 0)

    ######################

    # sum(parm$clust$C.m.vec) + parm$clust$C.m0 == parm$p



    ##########################################
    ##########################################
    ##### Computing proposal prob of reverse move
    ##########################################
    ##########################################

    old.prop <- 0

    for (gg in 1:old.parm$clust$G)
    {flag.gg <- (old.c.k==gg)
     count.gg <- sum(flag.gg)

     if (count.gg > 0)
     {old.prop <- old.prop + log(old.parm$clust$post.k[gg])*count.gg
     }
    }

    rho.prop <- new.prop - old.prop

    #######################################################
    #######################################################
    ########## computing true log-ratio: 2015
    #######################################################
    #######################################################

    # formula on page 1 of 07/27/15 notes

    # need to ensure there are no empty clusters in either
    # init.cc.parm or parm, otherwise likelihood formula
    # doesn't work (lgamma of negative values)

    tmp.new.parm <- parm
    indx.new <- parm$clust$C.m.vec > 0
    tmp.new.parm$clust$G <- sum(indx.new)
    tmp.new.parm$clust$C.m.vec <- parm$clust$C.m.vec[indx.new]
    #
    tmp.old.parm <- init.cc.parm
    indx.old <- init.cc.parm$clust$C.m.vec > 0
    tmp.old.parm$clust$G <- sum(indx.old)
    tmp.old.parm$clust$C.m.vec <- init.cc.parm$clust$C.m.vec[indx.old]

    rho.tru <- fn.d(d=parm$d, tmp.new.parm) - fn.d(d=parm$d, tmp.old.parm)

    new.log.lik <- 0


    for (gg in new.c.k)
    {indx.gg <- new.c.k==gg
		x_gg.mt <- matrix(parm$X[,I.k[indx.gg]], ncol=sum(indx.gg))
		##new.log.lik <- new.log.lik + sum(PDP_fn.log.lik(gg, x.mt=x_gg.mt, parm))
		new.log.lik <- new.log.lik + sum(PDP_fn.log.lik(gg, I.k[indx.gg], parm, data))
    }

    old.log.lik <- 0
    for (gg in old.c.k)
    {indx.gg <- old.c.k==gg
     x_gg.mt <- matrix(parm$X[,I.k[indx.gg]], ncol=sum(indx.gg))
     # old.log.lik <- old.log.lik + sum(PDP_fn.log.lik(gg, x.mt=x_gg.mt, old.parm))
	 old.log.lik <- old.log.lik + sum(PDP_fn.log.lik(gg, I.k[indx.gg], old.parm, data))
    }

    rho.tru <- rho.tru + new.log.lik - old.log.lik

    ########## toss a coin #################
	#print(rho.tru)
	#print(rho.prop)
    #print(new.log.lik)
	#print(old.log.lik)
	prob <- exp(min((rho.tru-rho.prop),0))
	#print(rho.tru)
	#print(rho.prop)
	
    flip<-as.logical(rbinom(n=1,size=1,prob=prob))

if (!flip) {parm <- init.cc.parm}


  } # end BIG if (!exit) loop


   list(parm, new.flag, exit, flip)
}


#####################################


PDP_fn.drop <- function(parm, computeMode)
 	{

	##########################################
	## Drop empty clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$G equal to number of non-empty clusters
	## (ii) Retain only clusters  1,...,parm$clust$G
	#########################################

	parm$clust$G <- sum(parm$clust$C.m.vec>0)
	num.dropped <- sum(parm$clust$C.m.vec==0)

  if (parm$clust$G > 0)
  {
	if (num.dropped > 0)
	{
	for (rr in 1:num.dropped)
		{
		old.label <-  min(which(parm$clust$C.m.vec==0))
		new.label <- max(which(parm$clust$C.m.vec>0))
		stopp <-  max(which(parm$clust$C.m.vec>0)) == parm$clust$G
		if (stopp)
			{break
			}
		parm <- PDP_fn.swap.clusters(parm, g1 = new.label, g2 = old.label, computeMode)
 		}
	}

	##########

	keep <- 1:parm$clust$G
	parm <- PDP_fn.clip.clusters(parm, keep)

	###########

	# parm$clust$K does not change (possibly some empty elementwise clusters)

	parm$N <- parm$n2 * parm$clust$G
	parm$clust$s.v <- as.vector(parm$clust$s.mt)

	parm$clust$n0 <- sum(parm$clust$s.v==0)
	parm$clust$n.vec <- array(,parm$clust$K)
	for (ss in 1:parm$clust$K)
		{parm$clust$n.vec[ss] <- sum(parm$clust$s.v==ss)  # HOT
		}

   }
   
  
	parm
	}


PDP_fn.orientation <- function(parm, cc_subset)
  {
  X.mt <- matrix(parm$X[,cc_subset], ncol=length(cc_subset))
  c.v <- parm$clust$c.v[cc_subset]
  orient.v <- parm$clust$orient.v[cc_subset]

  # manipulate PDP_fn.log.lik to get log-likelihoods separately for each sign
  tmp.parm <- parm
  tmp.parm$flip.sign=FALSE

  log_lik.mt <- array(, c(2, length(cc_subset)))

  for (gg in 0:parm$clust$G)
  {indx.gg <- which(c.v==gg)

    if (length(indx.gg)>0)
      {X_gg.mt <- matrix(X.mt[,indx.gg], ncol=length(indx.gg))
      log_lik.mt[1,indx.gg] <- PDP_fn.log.lik(gg, x.mt = X_gg.mt, tmp.parm)
      log_lik.mt[2,indx.gg] <- PDP_fn.log.lik(gg, x.mt = -X_gg.mt, tmp.parm)
      }
  }

  maxx.v <- apply(log_lik.mt, 2, max)

  log_lik.mt <- t(t(log_lik.mt) - maxx.v)
  lik.mt <- exp(log_lik.mt)


  # assuming equal prior prob to each sign

  for (tt in 1:length(cc_subset))
    {orient.v[tt] <- sample(c(-1,1),size=1,prob=lik.mt[,tt])
    }

  parm$clust$orient.v[cc_subset] <- orient.v

    parm
  }


fast_PDP_fn.main <- function(parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size, computeMode)
{
  p <- parm$p

  if (parm$standardize.X)
    {parm <- fn.standardize_orient.X(parm)
    }

	##########################
	# compute delta-neighborhoods
	#########################

  # SG: prob.compute.col.nbhd is the probability of computing
  # the neighborhoods, and col.frac.probes is the fraction of neighborhoods updated

  col_flip <- as.logical(rbinom(n=1,size=1,prob=prob.compute.col.nbhd))
  if (is.null(parm$clust$col.nbhd.k) | col_flip){
    parm <- PDP_fn.post.prob.and.delta(parm, data, max.col.nbhd.size, col.frac.probes, computeMode)
  }


	if (col.frac.probes < 1)
	{num_nbhds <- max(1,round(col.frac.probes*length(parm$clust$col.nbhd.k)))
	  parm$subset_nbhd.indx <- sort(sample(1:length(parm$clust$col.nbhd.k), size=num_nbhds))
	}

	if (col.frac.probes == 1)
	{parm$subset_nbhd.indx <- 1:length(parm$clust$col.nbhd.k)
	}

	new.flag.v <- NULL
    col.mh.flip.v <- col.mh.exit.v <- NULL

    for (cc in parm$subset_nbhd.indx)
	{previous.parm <- parm
	 parm$k <- parm$clust$col.nbhd.k[[cc]]

        if(length(parm$clust$col.nbhd[[cc]])==1)
		{ tmp <- PDP_fn.gibbs(k=parm$k, parm, data, computeMode)
		  parm <- tmp[[1]]
          new.flag.v <- c(new.flag.v, tmp[[2]])
		}

			if (length(parm$clust$col.nbhd[[cc]])>1)
			{
			  tmp <- PDP_fn.fast_col(cc, parm, data, computeMode)
			  parm <- tmp[[1]]
			  new.flag.v <- c(new.flag.v, tmp[[2]])
			  col.mh.exit.v <- c(col.mh.exit.v, tmp[[3]])
			  col.mh.flip.v <- c(col.mh.flip.v, tmp[[4]])
			}

	} # end for loop

	err <- PDP_fn.consistency.check(parm)
	if (err > 0)
			{stop(paste("LOOP: failed consistency check: err=",err))
	}

	parm$clust$col.new.flag <- mean(new.flag.v)
    if(!is.null(col.mh.flip.v))
    {parm$clust$col.mh.flip <- mean(col.mh.flip.v)
    parm$clust$col.mh.exit <- mean(col.mh.exit.v)
    } else{parm$clust$col.mh.flip <- 1
    parm$clust$col.mh.exit <- 1
    }

	##########################################
	## Drop empty group clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$G equal to number of non-empty clusters
	## (ii) Retain only clusters  1,...,parm$clust$G
	#########################################

	parm <- PDP_fn.drop(parm, computeMode)
	parm <- element_fn.drop(parm)

	err <- PDP_fn.compact.consistency.check(parm)
	if (err > 0)
		{stop(paste("MAIN FUNCTION END: failed consistency check: err=",err))
		}

	parm
	}



