
context("farma_acf")

ntest <- 60
test_that("farma_acf formula is correct.", {
  replicate(ntest, {
    # parameters
    N <- sample(10:20, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 2)
    nlag <- sample(1:3, 1)
    rho <- runif(nlag, -1, 1)
    # simplified formula
    acf1 <- farma_acf(alpha, 0, rho, dT, N)
    # long formula
    Tz <- Toeplitz(acf = fbm_acf(alpha, dT, N+nlag))
    J <- diag(c(rep(1 - sum(rho), N+nlag)))
    for(jj in 1:nlag) {
      J[cbind((jj+1):(N+nlag), 1:(N+nlag-jj))] <- rho[jj]  
    }
    T2 <- crossprod(J, Tz %*% J)[1:N,1:N]
    # check calculation
    max(abs(toeplitz(acf1) - T2))
    expect_equal(max(abs(toeplitz(acf1) - T2)), 0)
  })
})
