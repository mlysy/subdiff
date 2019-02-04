
context("fds_fit")

source("fit-functions.R")

# fbm downsampling loglikelihood
loglik <- function(theta, dX, dT, ds, Tz) {
  Xt <- dsample(apply(dX, 2, cumsum), ds)
  dX <- apply(Xt, 2, diff)
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- itrans_alpha(theta[1])
  mu <- theta[1+1:nd]
  Sigma <- itrans_Sigma(theta[1+nd+1:nq]) # default: log(D)
  Tz$setAcf(fbm_acf(alpha, dT*ds, N))
  suff <- lmn.suff(Y = dX, X = dT*ds, acf = Tz)
  lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
}

ntest <- 10
test_that("MLE is at the mode of the projection plots.", {
  skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
              "Package \"optimCheck\" required to run this test.")
  require(optimCheck)
  replicate(n = ntest, {
    # simulate data
    N <- sample(1000:2000, 1)
    dT <- runif(1)
    ds <- sample(2:10, 1)
    alpha <- runif(1, 0, 2)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rSnorm(n = ndims,
                           acf = fbm_acf(alpha, dT, N)))
    theta_hat <- fds_fit(dX, dT, ds, var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = floor(N/ds)-1) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik(theta, dX, dT, ds, Tz),
                         plot = FALSE, xrng = .05, npts = 20)
    expect_lt(max.xdiff(ocheck), .01)
  })
})
