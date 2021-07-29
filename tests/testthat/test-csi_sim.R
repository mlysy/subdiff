context("csi_sim")

test_that("csi_ZX is doing the same thing as the explicit calculation.", {
  cases <- expand.grid(ndim = 1:3,
                       nobs = c(2, 5, 10),
                       nsim = c(1, 3, 7))
  ntest <- nrow(cases)
  for(ii in 1:ntest) {
    ndim <- cases$ndim[ii]
    nobs <- cases$nobs[ii]
    nsim <- cases$nsim[ii]
    Sigma <- crossprod(matrix(rnorm(ndim^2),ndim,ndim))
    drift <- matrix(rnorm((nobs-1)*ndim), nobs-1, ndim)
    ## Sigma <- diag(rexp(ndim))
    ## drift <- matrix(0, nobs-1, ndim)
    Z <- matrix(rnorm((nobs-1)*ndim*nsim), nobs-1, ndim*nsim)
    X0 <- rnorm(ndim)
    Xt <- csi_ZX(Z,
                 drift = drift, Sigma = Sigma, X0 = X0,
                 nsim = nsim)
    Z2 <- array(Z, dim = c(nobs-1, nsim, ndim))
    Xt2 <- array(NA, dim = c(nobs, ndim, nsim))
    for(jj in 1:nsim) {
      dx <- matrix(Z2[,jj,], nobs-1, ndim) %*% chol(Sigma)
      Xt2[,,jj] <- apply(rbind(X0, dx + drift), 2, cumsum)
    }
    if(nsim == 1) Xt2 <- matrix(Xt2, nobs, ndim)
    expect_equal(Xt, Xt2)
  }
})
