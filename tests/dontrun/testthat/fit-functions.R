require(subdiff)
require(testthat)
require(optimCheck)

# dimension of covariacne matrix
getq <- function(ndims) if(ndims == 1) 1 else 3

# max of min of abs and rel error
max.xdiff <- function(x) {
  xdiff <- abs(diff(x))
  max(pmin(xdiff[,1], xdiff[,2]))
}

# different types of simulation
sim_func <- function(N, dT, type) {
  alpha <- runif(1, 0, 2)
  tau <- runif(1, 0, 1)
  rho <- runif(1, -1, 1)
  sigma2 <- runif(1, 0, 0.003)
  ndims <- sample(1:2, 1)
  switch(type, 
         "white" = matrix(rnorm(N*ndims), N, ndims), 
         "fBM" = as.matrix(rSnorm(n = ndims,
                                  acf = fbm_acf(alpha, dT, N))), 
         "fAR" = .far.sim(alpha, rho, dT, N, ndims),
         "fMA" = .fma.sim(alpha, rho, dT, N, ndims),
         "fdl" = .fdl.sim(alpha, tau, sigma2, dT, N, ndims))
}

.far.sim <- function(alpha, rho, dT, N, ndims) {
  dX <- as.matrix(rSnorm(n = ndims,
                         acf = fbm_acf(alpha, dT, N)))
  dY <- apply((1-rho)*dX, 2, function(x) {
    as.numeric(filter(x, filter = rho, method = "recursive"))
  })
  dY
}

.fma.sim <- function(alpha, rho, dT, N, ndims) {
  dX <- as.matrix(rSnorm(n = ndims,
                         acf = fbm_acf(alpha, dT, N+1)))
  dY <- (1-rho) * dX[-1,,drop=FALSE] + rho * dX[1:N,,drop=FALSE]
  dY
}

.fdl.sim <- function(alpha, tau, sigma2, dT, N, ndims) {
  dX <- as.matrix(rSnorm(n = ndims,
                         acf = fdyn_acf(alpha, tau, dT, N)))
  dZ <- matrix(rnorm((N+1)*ndims), N+1, ndims)
  dZ <- apply(dZ, 2, diff)
  dY <- dX + sqrt(sigma2) * dZ
  dY
}

