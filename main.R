#wkdir <- "C:\\Users\\Certified Copy\\Workspace\\Guha\\Laplace Approximation\\Variscan-Mixed-LXXY-OD"
#wkdir <- "C:\\Users\\cjfff_000\\Desktop\\Variscan-Mixed-LXXY-OD-18"
wkdir <- "C:\\Users\\Certified Copy\\Workspace\\Guha\\Laplace Approximation\\Variscan-Mixed-LXXY-OD-18"
setwd(wkdir)

library(msm)
library(cluster)

#' @import MASS
#' @import mvtnorm
#' @export
SimulateExample <- function(n = 25, p = 250, prop.data.type =c(0.5,0.5,0,0,0) , prop.X.miss=0, tau = 0.4419511, tau_0 = 1.25) {

	###################
	# generate covariates adding random noise of specified level
	# create objects data and true
	###################

	true_parm 			<- gen.clust(n, p)

	true_parm$tau 		<- tau
	true_parm$tau_0 	<- tau_0

	sim.X <- gen.X(n, p, prop.data.type, prop.X.miss, true_parm)

	simulation <- list(X = sim.X, parm = true_parm)
	class(simulation) <- "NPClustSimulation"

	return(simulation)
}


#' Fit an example DPP model
#'
#' @description \code{fitExample} fits an example DPP model
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' simulation <- simulateExample(n = 25, p = 125)
#'
#' # Fit model
#' posterior <- fitExample(simulation, n.burn = 10, n.reps = 20)
#'
#' # Summarize posterior
#' d_credible.v <- quantile(posterior$d.v, prob=c(.025,.975))
#' mean.taxicab <- mean(posterior$mean.taxicab.v)
#' se_mean.taxicab <- sd(posterior$mean.taxicab.v)/sqrt(length(posterior$mean.taxicab.v))
#' }
#'
#' @useDynLib NPCluster, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @export
fitExample <- function(data,
				n.burn = 10,
				n.reps = 20,
				max.row.nbhd.size = round(.1*25*125^.5), # should be small compared to n2*p^d (~ n2*G if d=.5)
				max.col.nbhd.size = round(.05*125), # should be small compared to p
				row.frac.probes = 0.05,
				col.frac.probes = .1,
                       		prob.compute.col.nbhd=.2,dahl.flag=FALSE,
				standardize.X=FALSE,
				flip.sign=FALSE, tBB_flag=FALSE,
				computeMode = createComputeMode()) {

	if (!inherits(data, "NPClustSimulation")) {
		stop("Wrong data structure")
	}

  if (!inherits(computeMode, "computeMode")) {
    stop("Wrong compute mode")
  }

  if (!standardize.X & flip.sign) {
    stop("Invalid input parameters-- flip.sign cannot be TRUE when standardize.X is FALSE")
  }


	###################
	# Detect clusters


	posterior <- fn.mcmc(text="CLUST ANALYZE...",
											 data$X$true, data$X$data,
											 n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes,
											 prob.compute.col.nbhd, data$parm, dahl.flag=dahl.flag, standardize.X, flip.sign, tBB_flag, computeMode)
	return (posterior)
}

#' @export
profileExample <- function(n = 25,
													 p = 250,
													 n.burn = 10,
													 n.reps = 20,
													 row.frac.probes = 0.05,
													 col.frac.probes = 0.05,
													 computeMode = createComputeMode(),
													 filename = "Rprof.out") {
	simulation <- simulateExample(n, p)

	Rprof(filename = filename, line.profiling = TRUE, interval = 0.001)
	posterior <- fitExample(simulation, n.burn = n.burn, n.reps = n.reps,
	           row.frac.probes = row.frac.probes,
	           col.frac.probes = col.frac.probes,
	           computeMode = computeMode)
	Rprof(NULL)
	#summaryRprof(lines = "show")$by.self
	return(posterior)
}


#' createComputeMode
#' @export
createComputeMode <- function(language = "R",
                              exactBitStream = FALSE,
                              extraSort = TRUE,
                              completeTest = FALSE,
                              tolerance = 1E-10,
                              test1 = FALSE,
                              test2 = FALSE,
                              test3 = FALSE) {
  if (!(language %in% c("C","R"))) {
    stop("Invalid language")
  }

  useR <- (language == "R")
  device <- NULL
  if (!useR) {
    doSort <- (exactBitStream | extraSort)
    device <- .createEngine(doSort)
  }

  object <- list(
    computeR = (language == "R" | completeTest),
    computeC = (language == "C"),
    device = device,
    exactBitStream = exactBitStream,
    extraSort = extraSort,
    tolerance = tolerance,
    test1 = test1,
    test2 = test2,
    test3 = test3
  )
  class(object) <- "computeMode"
  return(object)
}

#' assertEqual
assertEqual <- function(x, y, tolerance = 0) {
  if (length(x) != length(y)) {
    stop(cat("C++ error -- length:", length(x), length(y)))
  }
  if (any(abs(x - y) > tolerance)) {
    stop(cat("C++ error -- value:", x, y, tolerance, sep = "\n"))
  }
}


source("newDahl.R")
source("elementwise_DP.functions.R")
source("elementwise_main.R")
source("fast_PDP.functions.R")
source("gen.clust.R")
source("gen.X.R")
source("iterations.R")
source("lso.R")
source("NPCluster.R")
source("profile_code.R")
#source("profitable.R")
source("RcppExports.R")
source("split.merge.R")
source("variable.PDP.functions.R")



# Initial Paramters
n.burn = 500
n.reps = 1000
max.row.nbhd.size = round(.1*25*125^.5)
max.col.nbhd.size = round(.05*125)
row.frac.probes = 0.05
col.frac.probes = .1
prob.compute.col.nbhd=.2
dahl.flag=TRUE
standardize.X=FALSE
flip.sign=FALSE
tBB_flag=FALSE
computeMode = createComputeMode()

GBM <- read.table("Data\\GBM_Y2.csv",sep=",",header=TRUE)
GBM <- GBM[,-1]
data <- SimulateExample(n= 71,p = 393, prop.data.type =c(.07,0.365,0,0.165,0.40))
data$X$data$data.type <- c( rep(2,104), rep(4,86), rep(5,170), rep(1,33))
data$X$data$X  <- GBM
data$X$data$OX <- GBM


Run.m <-fitExample(data,
			n.burn = 5000,
			n.reps = 5000,
			max.row.nbhd.size = round(.1*25*125^.5), # should be small compared to n2*p^d (~ n2*G if d=.5)
			max.col.nbhd.size = round(.05*125), # should be small compared to p
			row.frac.probes = 0.05,
			col.frac.probes = .1,
                        prob.compute.col.nbhd=.2,
			dahl.flag=TRUE,
			standardize.X=FALSE,
			flip.sign=FALSE,
		        tBB_flag=FALSE,
			computeMode = createComputeMode())

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


