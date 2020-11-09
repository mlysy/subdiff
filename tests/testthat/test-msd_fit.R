
context("msd_fit")

test_that("msd calculation in C++ matches that of R", {
  test_cases <- expand.grid(ndims = 1:2,
                            demean = c(TRUE, FALSE))
  ntest <- nrow(test_cases)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(1000:2000, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 2)
    ndims <- test_cases$ndims[ii]
    dX <- as.matrix(rnormtz(n = ndims,
                            acf = fbm_acf(alpha, dT, N-1)))
    Xt <- apply(rbind(0, dX), 2, cumsum)
    # calculate msd
    nlag <- sample(nrow(Xt), 1)
    demean <- test_cases$demean[ii]
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
  }
})
