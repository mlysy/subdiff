
context("fma_acf")

ntest <- 20
test_that("fma_acf formula is correct.", {
  replicate(ntest, {
    # parameters
    N <- sample(10:20, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 2)
    rho <- runif(1, -1, 1)
    # simplified formula
    acf1 <- fma_acf(alpha, rho, dT, N)
    # long formula
    Tz <- Toeplitz(acf = fbm_acf(alpha, dT, N+1))
    J <- diag(c(rep(1 - rho, N+1)))
    J[cbind(1+1:N, 1:N)] <- rho
    T2 <- crossprod(J, Tz %*% J)[1:N,1:N]
    # check calculation
    expect_equal(max(abs(toeplitz(acf1) - T2)), 0)
  })
})
