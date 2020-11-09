
context("farma_acf")
source("fit-functions.R")

## ntest <- 60
test_that("farma_acf formula is correct for MA(q) filter.", {
  test_cases <- expand.grid(nlag = 0:4)
  ntest <- nrow(test_cases)
  for(ii in 1:ntest) {
    # parameters
    N <- sample(10:20, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 2)
    ## nlag <- sample(1:3, 1)
    nlag <- test_cases$nlag[ii]
    rho <- runif(nlag, -1, 1)/2
    # simplified formula
    acf1 <- farma_acf(alpha, rho = rho, dT = dT, N = N)
    # long formula
    acf2 <- fbm_acf(alpha, dT, N+nlag)
    if(nlag > 0) {
      Tz <- Toeplitz$new(acf = acf2)
      J <- diag(c(rep(1 - sum(rho), N+nlag)))
      for(jj in 1:nlag) {
        J[cbind((jj+1):(N+nlag), 1:(N+nlag-jj))] <- rho[jj]
      }
      T2 <- crossprod(J, Tz %*% J)[1:N,1:N]
      acf2 <- T2[1,]
    }
    # check calculation
    ## max(abs(toeplitz(acf1) - T2))
    expect_equal(max(abs(acf1 - acf2)), 0)
  }
})

test_that("farma_acf formula is correct for ARMA(1, q) filter.", {
  test_cases <- expand.grid(nlag = 0:4)
  ntest <- nrow(test_cases)
  for(ii in 1:ntest) {
    # parameters
    N <- sample(10:20, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 2)
    ## nlag <- sample(1:3, 1)
    nlag <- test_cases$nlag[ii]
    phi <- runif(1, -.5, .5)
    rho <- runif(nlag, -1, 1)/2
    # new formula
    acf1 <- farma_acf(alpha, phi = phi, rho = rho, dT = dT, N = N, m = 100)
    # old formula
    acf2 <- farma_acf2(alpha, phi = phi, rho = rho, dT, N = N, m = 100)
    # check calculation
    ## max(abs(toeplitz(acf1) - T2))
    expect_equal(max(abs(acf1 - acf2)), 0)
  }
})
