#' @keywords internal
gen.clust <- function(n, p) {

	true_parm <- NULL
	
	true_parm$d <- 0
	
	true_parm$clust <- NULL
	true_parm$clust$c.v <- array(,p)
	true_parm$clust$c.v[1] <- 1
	true_parm$clust$C.m.vec <- 1
	true_parm$clust$G <- 1
	true_parm$b0 <- 2.2
	true_parm$b1 <- 30
	
	
	for (xx in 2:p)
	{prob.v <- true_parm$clust$C.m.vec - true_parm$d
	 prob.v <- c(prob.v, (true_parm$b1 + true_parm$clust$G*true_parm$d))
	 new.c <- sample(1:(true_parm$clust$G+1), size=1, prob=prob.v)
	 
	 true_parm$clust$c.v[xx] <- new.c
	 new.flag <- (new.c > true_parm$clust$G)
	 
	 if (new.flag)
	 {true_parm$clust$G <- true_parm$clust$G + 1
	  true_parm$clust$C.m.vec <- c(true_parm$clust$C.m.vec, 1)
	 }
	 if (!new.flag)
	 {true_parm$clust$C.m.vec[new.c] <- true_parm$clust$C.m.vec[new.c] + 1
	 }
	}
	
	
	
	## neighborhood taxicab distances for column clusters
	
	tmp.mat <- array(0,c(p,p))
	
	for (jj in 1:true_parm$clust$G)
	{indx.jj <- which(true_parm$clust$c.v==jj)
	 tmp.mat[indx.jj,indx.jj] <- 1
	}
	
	true_parm$clust$nbhd.matrix <- tmp.mat
	
	#####
	
	true_parm$N <- n*true_parm$clust$G
	true_parm$clust$s.v <- array(,n*true_parm$clust$G)
	true_parm$clust$s.v[1] <- 1
	true_parm$clust$n.vec <- 1
	true_parm$clust$K <- 1
	true_parm$clust$M <- 11
	
	
	for (xx in 2:true_parm$N)
	{
	 prob.v <- c(true_parm$clust$n.vec)
	 prob.v <- c(prob.v, true_parm$clust$M)
	 new.s <- sample(1:(true_parm$clust$K+1), size=1, prob=prob.v)
	 
	 true_parm$clust$s.v[xx] <- new.s
	 new.flag <- (new.s > true_parm$clust$K)
	 
	 if (new.flag)
	 {true_parm$clust$K <- true_parm$clust$K + 1
	  true_parm$clust$n.vec <- c(true_parm$clust$n.vec, 1)
	 }

	 if ((!new.flag))
	 {true_parm$clust$n.vec[new.s] <- true_parm$clust$n.vec[new.s] + 1
	 }
	}
	
	
	true_parm$clust$s.mt <- matrix(true_parm$clust$s.v, nrow=n)
	
	###########
	
	true_parm$clust$mu2  <- 0.193035
	true_parm$clust$tau2 <- 1.179075
	
	true_parm$clust$phi.v <- rnorm(n=true_parm$clust$K, mean=true_parm$clust$mu2, sd=true_parm$clust$tau2)
	

	
	
	
	return (true_parm)
}
