#Old Results
#Kmediod Results
#n=85,p=300,cont.5, binary d.25, binary ind.25 |||.03364  ||| 0.07
#n=85,p=300,cont.5, binary true .5 |||.02331  ||| 0.014
#n=85,p=300,cont.5, binary dep .5 |||.0132  ||| 0.013 (nearly the same)

#Comparison Code
#Scenario cases
#d=0,0.2,0.4,0.6,0.8

data		 <- SimulateExample(n= 85,p = 299, prop.data.type =c(0.1,0.5,0,0.2,0.2))

#K-mediods
library(cluster)

dist <- daisy(t(data$X$data$X),metric="gower")

#Select K

sil_width <- c(NA)


for(i in 1:49 ){
 pam_fit 	  <- pam(dist,diss = TRUE, k = 5*i)
 sil_width[i] <- pam_fit$silinfo$avg.width
}

plot(1:49, sil_width, xlab="Number of Clusters",ylab ="Silhoutte Width")
lines(1:49,sil_width)

pam_fit <- pam(dist, diss = TRUE, k =5*49)
mat.pam <- matrix(0,299,299)

for(i in 1:299){
 mat.pam[which(pam_fit$clustering == pam_fit$clustering[i] ),i ] <- 1
}

mean(mat.pam != data$parm$clust$nbhd.matrix)

