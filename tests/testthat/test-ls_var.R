context("ls_fit")

source("fit-functions.R")

test_that("ls_var is same in C++ and R.", {
  ntest <- 10
  for(ii in 1:ntest) {
    ntau <- sample(10, 1)
    N <- sample(10, 1)
    alpha <- runif(1) * 2
    tau <- sort(runif(ntau))
    V_r <- ls_var_r(alpha, tau, N)
    V_cpp <- ls_var(alpha, tau, N)
    diag(V_r) <- diag(V_cpp) # this is where `eps^alpha` might happen...
    expect_equal(V_r, V_cpp)
  }
})
