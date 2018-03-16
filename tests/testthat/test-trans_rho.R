
context("rho trans/itrans")

ntest <- 20

test_that("rho = itrans(trans(rho))", {
  for(ii in 1:ntest) {
    rho <- runif(1,-1,1)
    nu <- trans_rho(rho)
    expect_equal(rho, itrans_rho(nu))
  }
})

test_that("nu = trans(itrans(nu))", {
  for(ii in 1:ntest) {
    nu <- rnorm(1)
    rho <- itrans_rho(nu)
    expect_equal(nu, trans_rho(rho))
  }
})
