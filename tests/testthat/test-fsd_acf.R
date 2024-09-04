context("fsd_acf")

## source("subdiff-testfunctions.R")

ntest <- 10
test_that("fsd_acf formula is correct.", {
  replicate(ntest, {
    # parameters
    N <- sample(10:20, 1)
    dt <- runif(1)
    tau <- 1/sample(5:10, size = 1)
    M <- 50 # resolution
    alpha <- runif(1, 0, 2)
    # simplified formula
    acf1 <- fsd_acf(alpha, tau, 0, dt, N)
    # long formula
    acf_2 <- fbm_acf(alpha, dt*tau, N*M/tau)
    trans_mat <- function(N, M, tau) {
      mat <- matrix(0, N, N*M/tau)
      for(ii in 1:N){
        mat[ii, (ii-1)*M/tau + 1:M] <- 1/M
      }
      mat
    }
    acf2 <- trans_mat(N, M, tau) %*% as.matrix(acf_2)
    # check calculation
    expect_lt(max_xdiff(cbind(acf1, acf2)), .01)
  })
})


