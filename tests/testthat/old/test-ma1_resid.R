
context("ma1_resid")

ntest <- 20
test_that("Z = resid(filter(Z))", {
  for(jj in 1:ntest) {
    N <- sample(10:50, 1)
    qq <- sample(1:5, 1)
    rho <- runif(1,-.5,.5)
    Z <- matrix(rnorm(N*qq), N, qq)
    X <- matrix(NA, N, qq)
    X[1,] <- (1-rho) * Z[1,]
    for(ii in 2:N) {
      X[ii,] <- (1-rho) * Z[ii,] + rho * Z[ii-1,]
    }
    Z2 <- ma1_resid(as.matrix(X), rho)
    expect_equal(Z, Z2, tol = 1e-5)
  }
})
