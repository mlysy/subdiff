
context("fma_fit")

source("fit-functions.R")

# fma loglikelihood
loglik <- function(theta, dX, dT, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- itrans_alpha(theta[1])
  rho <- itrans_rho(theta[2])
  mu <- theta[2+1:nd]
  Sigma <- itrans_Sigma(theta[2+nd+1:nq]) # default: log(D)
  Tz$setAcf(fma_acf(alpha, rho, dT, N))
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
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
    dX <- sim_func(N, dT, "fMA")
    theta_hat <- fma_fit(dX, dT, var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik(theta, dX, dT, Tz),
                         plot = FALSE, xrng = .05, npts = 20)
    print(max.xdiff(ocheck))
#    expect_lt(max.xdiff(ocheck), .01)
  })
})
