require(subdiff)

context("alpha trans/itrans")

ntest <- 20

test_that("alpha = itrans(trans(alpha))", {
  for(ii in 1:ntest) {
    alpha <- 2*runif(1)
    gamma <- trans_alpha(alpha)
    expect_equal(alpha, itrans_alpha(gamma))
  }
})

test_that("gamma = trans(itrans(gamma))", {
  for(ii in 1:ntest) {
    gamma <- rnorm(1)
    alpha <- itrans_alpha(gamma)
    expect_equal(gamma, trans_alpha(alpha))
  }
})
