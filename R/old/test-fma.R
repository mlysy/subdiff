require(subdiff)
require(doParallel) # parallelize


# functions ---------------------------------------------------------------

ma_seq <- function(dX, rho) {
  N <- nrow(dX)
  d <- ncol(dX)
  if(length(rho) == 1) {
    dY <- (1-rho)*dX + rho*rbind(0, dX[-N,])
  } else if(length(rho) == 2) {
    dY <- (1-sum(rho))*dX + rho[1]*rbind(0, dX[-N,]) + rho[2]*rbind(matrix(0,2,d), dX[-c(N-1,N),])
  } else if(length(rho == 3)) {
    dY <- (1-sum(rho))*dX + rho[1]*rbind(0, dX[-N,]) + rho[2]*rbind(matrix(0,2,d), dX[-c(N-1,N),]) + 
      rho[3]*rbind(matrix(0,3,d), dX[-c(N-2,N-1,N),])
  }
}

coverage_func <- function(data_stats, value, ncoef) {
  nrep <- length(data_stats)
  ncount <- rep(NA, nrep)
  for(ii in 1:nrep) {
    stat <- data_stats[[ii]]
    if(!is.list(stat)) {
      ncount[ii] <- 0
    } else {
      mu <- stat$coef[ncoef]
      sigma <- sqrt(stat$vcov[ncoef,ncoef])
      if(mu-1.96*sigma < value && mu+1.96*sigma > value) {
        ncount[ii] <- 1
      } else {
        ncount[ii] <- 0
      } 
    }
  }
  ncount
}


# simulation and inference ------------------------------------------------

rho <- c(.2, .1, -.5)
acf1 <- fbm_acf(.8, 1/60, 2e3)
nrep <- 1e3

nclust <- detectCores(logical = FALSE)
cl <- makeCluster(nclust)
registerDoParallel(cl)
clusterSetRNGStream(cl = cl)
cl_pkg <- c("subdiff", "nortest")

stats <- foreach(ii = 1:nrep, .packages = cl_pkg) %dopar% {
  dX <- rSnorm(n = 2, acf = acf1)
  dY <- ma_seq(dX, rho)
  fit <- try(fma_fit(dX = dY, dT = 1/60, nlag = 3, var_calc = TRUE), silent = TRUE)
  fit
}

stopCluster(cl)

sum(coverage_func(stats, trans_alpha(.8), "gamma")) / nrep


