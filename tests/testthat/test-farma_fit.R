
context("farma_fit")

source("fit-functions.R")

# farma(1,1) loglikelihood
loglik_11 <- function(theta, dX, dt, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- ilogit(theta[1], min = 0, max = 2)
  phi <- ilogit(theta[2], min = -1, max = 1)
  rho <- ilogit(theta[3], min = -1, max = 1)
  mu <- theta[3+1:nd]
  Sigma <- itrans_Sigma(theta[3+nd+1:nq]) # default: log(D)
  Tz$set_acf(farma_acf(alpha, phi, rho, dt, N))
  suff <- lmn_suff(Y = dX, X = dt, V = Tz, Vtype = "acf")
  lmn_loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
}

ntest <- 1
test_that("MLE is at the mode of the projection plots.", {
  skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
              "Package \"optimCheck\" required to run this test.")
  require(optimCheck)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(1000:2000, 1)
    dt <- runif(1)
    alpha <- runif(1, 0, 2)
    rho <- runif(1, -1, 1)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rnormtz(n = ndims,
                            acf = fbm_acf(alpha, dt, N+1)))
    dX <- (1-rho) * dX[-1,,drop=FALSE] + rho * dX[1:N,,drop=FALSE]
    model <- farma_model$new(dX, dt, p = 1, q = 1)
    theta_hat <- farma_fit(dX, dt, order = c(1,1), vcov = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz$new(N = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik_11(theta, dX, dt, Tz),
                         plot = F, xrng = .1, npts = 20)
    expect_lt(max_xdiff(ocheck), .05)
  }
})

context("fma_fit")


# farma(0,1) loglikelihood
loglik_01 <- function(theta, dX, dt, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- ilogit(theta[1], min = 0, max = 2)
  rho <- ilogit(theta[2], min = -1, max = 1)
  mu <- theta[2+1:nd]
  Sigma <- itrans_Sigma(theta[2+nd+1:nq]) # default: log(D)
  Tz$set_acf(farma_acf(alpha, numeric(), rho, dt, N))
  suff <- lmn_suff(Y = dX, X = dt, V = Tz, Vtype = "acf")
  lmn_loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
}

ntest <- 1
test_that("MLE is at the mode of the projection plots.", {
  skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
              "Package \"optimCheck\" required to run this test.")
  require(optimCheck)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(1000:2000, 1)
    dt <- runif(1)
    alpha <- runif(1, 0, 2)
    rho <- runif(1, -1, 1)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rnormtz(n = ndims,
                            acf = fbm_acf(alpha, dt, N+1)))
    dX <- (1-rho) * dX[-1,,drop=FALSE] + rho * dX[1:N,,drop=FALSE]
    theta_hat <- farma_fit(dX, dt, order = c(0,1), vcov = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz$new(N = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik_01(theta, dX, dt, Tz),
                         plot = F, xrng = .1, npts = 20)
    expect_lt(max_xdiff(ocheck), .05)
  }
})

