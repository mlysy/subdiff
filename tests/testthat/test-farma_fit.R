
context("farma_fit")

source("fit-functions.R")

# fma loglikelihood
loglik <- function(theta, dX, dT, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- itrans_alpha(theta[1])
  gamma <- itrans_gamma(theta[2])
  rho <- itrans_rho(theta[3])
  mu <- theta[3+1:nd]
  Sigma <- itrans_Sigma(theta[3+nd+1:nq]) # default: log(D)
  Tz$setAcf(farma_acf(alpha, gamma, rho, dT, N))
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
    alpha <- runif(1, 0, 2)
    rho <- runif(1, -1, 1)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rSnorm(n = ndims,
                           acf = fbm_acf(alpha, dT, N+1)))
    dX <- (1-rho) * dX[-1,,drop=FALSE] + rho * dX[1:N,,drop=FALSE]
    theta_hat <- farma_fit(dX, dT, var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik(theta, dX, dT, Tz),
                         plot = FALSE, xrng = .1, npts = 20)
    expect_lt(max.xdiff(ocheck), .05)
  })
})


# fma loglikelihood
loglik1 <- function(theta, dX, dT, alpha, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  gamma <- itrans_gamma(theta[1])
  rho <- itrans_rho(theta[2])
  mu <- theta[2+1:nd]
  Sigma <- itrans_Sigma(theta[2+nd+1:nq]) # default: log(D)
  Tz$setAcf(farma_acf(alpha, gamma, rho, dT, N))
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
}


test_that("MLE is at the mode of the projection plots.", {
  skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
              "Package \"optimCheck\" required to run this test.")
  require(optimCheck)
  replicate(n = ntest, {
    # simulate data
    N <- sample(1000:2000, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 2)
    rho <- runif(1, -1, 1)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rSnorm(n = ndims,
                           acf = fbm_acf(alpha, dT, N+1)))
    dX <- (1-rho) * dX[-1,,drop=FALSE] + rho * dX[1:N,,drop=FALSE]
    theta_hat <- farma_fit(dX, dT, alpha = alpha, var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik1(theta, dX, dT, alpha, Tz),
                         plot = FALSE, xrng = .1, npts = 20)
    expect_lt(max.xdiff(ocheck), .05)
  })
})
