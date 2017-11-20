require(subdiff)

context("Sigma trans/itrans")

ntest <- 20
cases <- expand.grid(run = 1:ntest, q = 1:2,
                     D_scale = c(TRUE, FALSE), rho_scale = c(TRUE, FALSE))
ncases <- nrow(cases)

test_that("Sigma = itrans(trans(Sigma))", {
  for(ii in 1:ncases) {
    q <- cases$q[ii]
    D_scale <- cases$D_scale[ii]
    rho_scale <- cases$rho_scale[ii]
    if(rho_scale) {
      Sigma <- if(q == 1) rexp(1) else c(rexp(2), runif(1,-1,1))
    } else {
      Sigma <- if(q == 1) t(rexp(1)) else crossprod(matrix(rnorm(4),2,2))
    }
    lambda <- trans_Sigma(Sigma, D_scale = D_scale, rho_scale = rho_scale)
    expect_equal(Sigma, itrans_Sigma(lambda, D_scale = D_scale,
                                     rho_scale = rho_scale))
  }
})

test_that("lambda = trans(itrans(lambda))", {
  for(ii in 1:ncases) {
    q <- cases$q[ii]
    D_scale <- cases$D_scale[ii]
    rho_scale <- cases$rho_scale[ii]
    lambda <- if(q == 1) rnorm(1) else rnorm(3)
    Sigma <- itrans_Sigma(lambda, D_scale = D_scale, rho_scale = rho_scale)
    expect_equal(lambda, trans_Sigma(Sigma, D_scale = D_scale,
                                     rho_scale = rho_scale))
  }
})
