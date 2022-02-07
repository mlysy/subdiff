context("trans_Sigma and itrans_Sigma")

source("subdiff-testfunctions.R")

test_that("Sigma = itrans(trans(Sigma))", {
  ntest <- 20
  for(ii in 1:ntest) {
    ndim <- sample(1:10, 1)
    Sigma <- crossprod(matrix(rnorm(ndim^2), ndim, ndim))
    lambda <- trans_Sigma(Sigma)
    expect_equal(Sigma, itrans_Sigma(lambda))
  }
})

test_that("lambda = trans(itrans(lambda))", {
  ntest <- 20
  for(ii in 1:ntest) {
    ndim <- sample(1:10, 1)
    lambda <- rnorm(ndim*(ndim+1)/2)
    Sigma <- itrans_Sigma(lambda)
    expect_equal(lambda, trans_Sigma(Sigma))
  }
})
