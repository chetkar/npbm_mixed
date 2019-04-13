
fn.Dahl <- function(parm, All.Stuff, dahl.flag)
{### Dahl calculations
  
  tmp.mat <- array(0,c(parm$p,parm$p))
  
  for (jj in 1:parm$clust$G)
  {indx.jj <- which(parm$clust$c.v==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  parm$dahlDist.mt = tmp.mat
  
  if (dahl.flag)
  {# interpretation: square of this is double the total number of mismatching pairs of allocations
    parm$Dahl_dist =  sqrt(sum((parm$dahlDist.mt - All.Stuff$meanDahlDist.mt)^2))
  
  if (parm$Dahl_dist < parm$min_Dahl_dist)
  {parm$min_Dahl_dist = parm$Dahl_dist
  parm$est_c.v  <- parm$clust$c.v
  }
  }
  
  parm
}
