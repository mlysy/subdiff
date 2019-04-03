
context("fbm_fit")

source("fit-functions.R")

# fbm loglikelihood
loglik <- function(theta, dX, dT, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- ilogit(theta[1], min = 0, max = 2)
  mu <- theta[1+1:nd]
  Sigma <- itrans_Sigma(theta[1+nd+1:nq]) # default: log(D)
  Tz$setAcf(fbm_acf(alpha, dT, N))
  suff <- lmn.suff(Y = dX, X = dT, V = Tz, Vtype = "acf")
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
    alpha <- runif(1, 0, 2)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rSnorm(n = ndims,
                           acf = fbm_acf(alpha, dT, N)))
    theta_hat <- fbm_fit(dX, dT, var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik(theta, dX, dT, Tz),
                         plot = FALSE, xrng = .05, npts = 20)
    expect_lt(max.xdiff(ocheck), .01)
  })
})
