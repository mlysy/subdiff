context("log1pe")

source("subdiff-testfunctions.R")

test_that("log1pe is same in C++ and R.", {
  ntest <- 10
  for(ii in 1:ntest) {
    N <- sample(10, 1)
    A <- runif(1, 0, log(100))
    x <- runif(N, -exp(A), exp(A))
    expect_equal(subdiff:::log1pe(x), log1pexp(x), tolerance = 1e-6)
  }
})
