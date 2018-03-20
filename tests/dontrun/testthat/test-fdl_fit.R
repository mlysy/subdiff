
context("fbm_fit")

source("fit-functions.R")

# fdl loglikelihood
loglik <- function(theta, dX, dT, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- itrans_alpha(theta[1])
  tau <- itrans_tau(theta[2])
  sigma2 <- theta[3]^2
  mu <- theta[3+1:nd]
  Sigma <- itrans_Sigma(theta[3+nd+1:nq]) # default: log(D)
  Tz$setAcf(fdyn_acf(alpha, tau, dT, N) + sigma2 * c(2,-1,rep(0,N-2)))
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
    dX <- sim_func(N, dT, "white")
    theta_hat <- fdl_fit(dX, dT, var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik(theta, dX, dT, Tz),
                         plot = FALSE, xrng = .05, npts = 20)
    print(max.xdiff(ocheck))
    optim_proj(xsol = theta_hat,
               fun = function(theta) loglik(theta, dX, dT, Tz),
               plot = TRUE, xrng = .05, npts = 20)
#     expect_lt(max.xdiff(ocheck), .01)
  })
})
