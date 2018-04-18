
context("msd_fit")

ntest <- 12

test_that("msd calculation in C++ matches that of R", {
  replicate(ntest, {
    # simulate data
    N <- sample(1000:2000, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 2)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rSnorm(n = ndims,
                           acf = fbm_acf(alpha, dT, N-1)))
    Xt <- apply(rbind(0, dX), 2, cumsum)
    # calculate msd
    nlag <- sample(nrow(Xt), 1)
    demean <- rbinom(1,1,.5) == 1
    # C++ estimate
    msd_C <- msd_fit(Xt, nlag, demean)
    # R estimate
    if(demean) {
      mu <- colMeans(dX)
      Zt <- Xt - (1:N-1) %o% mu
    } else {
      Zt <- Xt
    }
    msd_R <- rep(NA, nlag)
    for(ii in 1:nlag) {
      msd_R[ii] <- mean(apply(Zt, 2, diff, lag = ii)^2)
    }
    expect_equal(max(abs(msd_C - msd_R)), 0)
  })
})
