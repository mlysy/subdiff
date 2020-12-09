context("ls_fit")

source("fit-functions.R")

test_that("ls_var is same in C++ and R.", {
  ntest <- 10
  for(ii in 1:ntest) {
    nlags <- sample(10, 1)
    N <- sample(100:200, 1)
    alpha <- runif(1, .25, 1.5)
    lags <- sort(sample(20, 5))
    V_r <- ls_var_r(alpha, lags, N)
    V_cpp <- subdiff:::ls_var(alpha, lags, N)
    diag(V_r) <- diag(V_cpp) # this is where `eps^alpha` might happen...
    expect_equal(V_r, V_cpp, tolerance = 1e-6)
  }
})
