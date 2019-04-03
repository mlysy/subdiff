
context("farma_fit")

source("fit-functions.R")

# farma(1,1) loglikelihood
loglik_11 <- function(theta, dX, dT, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- ilogit(theta[1], min = 0, max = 2)
  phi <- ilogit(theta[2], min = -1, max = 1)
  rho <- ilogit(theta[3], min = -1, max = 1)
  mu <- theta[3+1:nd]
  Sigma <- itrans_Sigma(theta[3+nd+1:nq]) # default: log(D)
  Tz$setAcf(farma_acf(alpha, phi, rho, dT, N))
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
    rho <- runif(1, -1, 1)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rSnorm(n = ndims,
                           acf = fbm_acf(alpha, dT, N+1)))
    dX <- (1-rho) * dX[-1,,drop=FALSE] + rho * dX[1:N,,drop=FALSE]
    theta_hat <- farma_fit(dX, dT, order = c(1,1), var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik_11(theta, dX, dT, Tz),
                         plot = F, xrng = .1, npts = 20)
    expect_lt(max.xdiff(ocheck), .05)
  })
})

# farma(0,1) loglikelihood
loglik_01 <- function(theta, dX, dT, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- ilogit(theta[1], min = 0, max = 2)
  rho <- ilogit(theta[2], min = -1, max = 1)
  mu <- theta[2+1:nd]
  Sigma <- itrans_Sigma(theta[2+nd+1:nq]) # default: log(D)
  Tz$setAcf(farma_acf(alpha, 0, rho, dT, N))
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
    rho <- runif(1, -1, 1)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rSnorm(n = ndims,
                           acf = fbm_acf(alpha, dT, N+1)))
    dX <- (1-rho) * dX[-1,,drop=FALSE] + rho * dX[1:N,,drop=FALSE]
    theta_hat <- farma_fit(dX, dT, order = c(0,1), var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik_01(theta, dX, dT, Tz),
                         plot = F, xrng = .1, npts = 20)
    expect_lt(max.xdiff(ocheck), .05)
  })
})

