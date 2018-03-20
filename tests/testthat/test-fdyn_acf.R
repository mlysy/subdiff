
context("fdyn_acf")

ntest <- 20
test_that("fdyn_acf formula is correct.", {
  replicate(ntest, {
    # parameters
    N <- sample(10:20, 1)*10
    nrep <- 1e2
    dT <- 1/60
    tau <- 1/10
    M <- 50
    alpha <- .8
    # simplified formula
    acf1 <- fdyn_acf(alpha, tau, dT, 30)
    # long formula
    dY <- rSnorm(nrep, acf = fbm_acf(alpha, dT*tau/M, N/tau*M))
    Yt <- apply(dY, 2, cumsum)
    Xt <- matrix(NA, N, nrep)
    for(ii in 1:nrep) {
      for(jj in 1:N) {
        Xt[jj, ii] <- mean(Yt[((jj-1)/tau*M) + 1:M,ii])
      }
    }
    dX <- apply(Xt, 2, diff)
    acf_mat <- sapply(1:nrep, function(ii) {
      acf(dX[,ii], plot = FALSE, type = "covariance", lag.max = 29)$acf
    })
    acf2 <- apply(acf_mat, 1, mean)
    # check calculation
    expect_lt(max.xdiff(cbind(acf1, acf2)), .05)
  })
})
