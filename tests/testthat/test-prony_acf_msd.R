
context("Prony-GLE ACF and MSD")

## # maximum absolute and relative error between two arrays
## max.diff <- function(x1, x2) {
##   c(abs = max(abs(x1-x2)), rel = max(abs(x1-x2)/max(abs(x1), 1e-8)))
## }

ntest <- 1
test_that("MSD => ACF", {
  for(ii in 1:ntest) {
    Temp <- runif(1, 290, 300)
    vsigma <- runif(1, 0.01, 0.02)
    K <- sample(10:200, 1)
    alpha <- runif(1, .1, .9)
    tau <- runif(1, 1e-4, 1e-2)
    lambda <- ((1:K)/K)^(1/alpha)/tau
    dt <- runif(1)
    N <- 1e4
    tseq <- 1:N*dt
    acf <- prony_acf(lambda, vsigma, dt, N)
    msd <- prony_msd(tseq, lambda, vsigma)
    expect_equal(acf, msd2acf(msd))
  }
})

test_that("ACF => MSD", {
  for(ii in 1:ntest) {
    Temp <- runif(1, 290, 300)
    vsigma <- runif(1, 0.01, 0.02)
    K <- sample(10:200, 1)
    alpha <- runif(1, .1, .9)
    tau <- runif(1, 1e-4, 1e-2)
    lambda <- ((1:K)/K)^(1/alpha)/tau
    dt <- runif(1)
    N <- 1e4
    tseq <- 1:N*dt
    acf <- prony_acf(lambda, vsigma, dt, N)
    msd <- prony_msd(tseq, lambda, vsigma)
    expect_equal(msd, acf2msd(acf))
  }
})
