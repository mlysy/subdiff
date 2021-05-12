#--- check the LS var formula --------------------------------------------------

require(subdiff)

alpha <- .5
dt <- runif(1, max = .1)
N <- 2000
ndim <- 3
nreps <- 1000

acf <- fbm_acf(alpha = alpha, dt = dt, N = N)
dX <- SuperGauss::rnormtz(n = ndim*nreps, acf = acf)
dX <- array(dX, dim = c(N, ndim, nreps))

msd <- apply(dX, 3, function(dx) {
  Xt <- apply(dx, 2, cumsum)
  msd_fit(Xt, demean = FALSE)
})
msd <- t(msd)

nlags <- ncol(msd)
V1 <- var(log(msd))
## V2 <- (alpha/ndim) * subdiff:::ls_var(alpha = 1, lags = 1:nlags, N = N)
V2 <- (1/ndim) * subdiff:::ls_var(alpha = alpha, lags = 1:nlags, N = N)

ind <- sample(nlags, 2, replace = TRUE)
## (V1[ind[1], ind[2]] - V2[ind[1], ind[2]]) / V1[ind[1], ind[2]]
message("[", ind[1], ",", ind[2], "] = ",
        signif(sqrt(V1[ind[1], ind[2]]/V2[ind[1], ind[2]]), 3))


#--- scratch -------------------------------------------------------------------

msd <- replicate(nreps, {
  Xt <- cumsum(SuperGauss::rnormtz(acf = acf))
  msd_fit(matrix(Xt), demean = FALSE)
})
